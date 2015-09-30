#!/usr/bin/env python2.7

# This file is part of pyGenClean.
#
# pyGenClean is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# pyGenClean is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# pyGenClean.  If not, see <http://www.gnu.org/licenses/>.


import os
import re
import sys
import logging
import argparse

import numpy as np

import matplotlib as mpl

from .. import __version__


logger = logging.getLogger("plot_eigenvalues")


def main(argString=None):
    """The main function.

    The purpose of this module is to plot Eigenvectors provided by the
    Eigensoft software.

    Here are the steps of this module:

    1. Reads the Eigenvector (:py:func:`read_eigenvalues`).
    2. Plots the Scree Plot (:py:func:`create_scree_plot`).

    """
    # Getting and checking the options
    args = parse_args(argString)
    check_args(args)

    # Reads the eigenvalues
    eigenvalues = read_eigenvalues(args.evec)

    # Creates the plot
    create_scree_plot(eigenvalues, args.out, args)


def create_scree_plot(data, o_filename, options):
    """Creates the scree plot.

    :param data: the eigenvalues.
    :param o_filename: the name of the output files.
    :param options: the options.

    :type data: numpy.ndarray
    :type o_filename: str
    :type options: argparse.Namespace

    """
    # Importing plt
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    plt.ioff()

    # Computing the cumulative sum
    cumul_data = np.cumsum(data)

    # Creating the figure and axes
    fig, axes = plt.subplots(2, 1, figsize=(8, 16))

    # The title of the figure
    fig.suptitle(options.scree_plot_title, fontsize=16, weight="bold")

    fig.subplots_adjust(hspace=0.27, top=0.93, bottom=0.06)

    # Modifying the spines
    for axe in axes:
        axe.xaxis.set_ticks_position("bottom")
        axe.yaxis.set_ticks_position("left")
        axe.spines["top"].set_visible(False)
        axe.spines["right"].set_visible(False)
        axe.spines["left"].set_position(("outward", 9))
        axe.spines["bottom"].set_position(("outward", 9))

    # First, plotting the eigenvalues
    axes[0].set_title("Scree Plot", weight="bold")
    axes[0].set_xlabel("Component")
    axes[0].set_ylabel("Eigenvalue")
    axes[0].plot(np.arange(len(data)) + 1, data, marker="o", c="#0099CC",
                 mec="#0099CC", ls="-", lw=2, clip_on=False)

    # Then, plotting the annotation
    for i in range(len(data)):
        axes[0].annotate(np.round(data[i], 3),
                         xy=(np.arange(len(data))[i] + 1, data[i]),
                         xytext=(1, 10), textcoords="offset points",
                         ha="left", va="bottom",
                         bbox=dict(boxstyle="round,pad=0.5", fc="#FFFFFF"))

    # Next plot the cumulative values
    axes[1].set_title(("Cumulative explained variance "
                       "(max={:.3f})".format(np.sum(data))), weight="bold")
    axes[1].set_xlabel("Component")
    axes[1].set_ylabel("Cumulative explained variance")
    axes[1].axhline(np.sum(data) * 0.8, ls="--", lw="2", c="#999999")
    axes[1].plot(np.arange(len(data)) + 1, cumul_data, marker="o", c="#CC0000",
                 mec="#CC0000", mfc="#CC0000", ls="-", lw=2, clip_on=False)

    # Then, plotting the annotation
    for i in range(len(data)):
        axes[1].annotate(np.round(cumul_data[i], 3),
                         xy=(np.arange(len(data))[i] + 1, cumul_data[i]),
                         xytext=(1, -10), textcoords="offset points",
                         ha="left", va="top",
                         bbox=dict(boxstyle="round,pad=0.5", fc="#FFFFFF"))

    # Saving the file
    plt.savefig(o_filename, dpi=300)
    plt.close(fig)


def read_eigenvalues(i_filename):
    """Reads the eigenvalues from EIGENSOFT results.

    :param i_filename: the name of the input file.

    :type i_filename: str

    :returns: a :py:class:`numpy.ndarray` array containing the eigenvalues.

    """
    # The data is the first line of the result file (should begin with
    # "#eigvals"
    data = None
    with open(i_filename, "r") as i_file:
        data = re.split(
            r"\s+",
            re.sub(r"(^\s+)|(\s+$)", "", i_file.readline()),
        )
        if not data[0].startswith("#eigvals"):
            m = "{}: not a evec file".format(i_filename)
            raise ProgramError(m)
        data = np.array(data[1:], dtype=float)

    return data


def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    # Checking that the input file exists
    if not os.path.isfile(args.evec):
        m = "{}: no such file".format(args.evec)
        raise ProgramError(m)
    return True


def parse_args(argString=None):
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ====================== ====== ================================
            Options         Type            Description
    ====================== ====== ================================
    ``--evec``             string The EVEC file from EIGENSOFT
    ``--scree-plot-title`` string The main title of the scree plot
    ``--out``              string The name of the output file
    ====================== ====== ================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    args = None
    if argString is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argString)

    return args


def add_custom_options(parser):
    """Adds custom options to a parser.

    :param parser: the parser to which the options will be added.

    :type parser: argparse.ArgumentParser

    """
    parser.add_argument("--scree-plot-title", type=str, metavar="TITLE",
                        default="EIGENSOFT results",
                        help="The main title of the scree plot [%(default)s]")


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.

    :param msg: the message to print to the user before exiting.

    :type msg: str

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user.

        :type msg: str

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


# The parser object
desc = "Plots eigenvalues"
long_desc = None
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The input files
group = parser.add_argument_group("Input File")
group.add_argument("--evec", type=str, metavar="FILE", required=True,
                   help="The EVEC file from EIGENSOFT")

# The plot options
group = parser.add_argument_group("Plot Options")
add_custom_options(group)

# The output
group = parser.add_argument_group("Output Options")
group.add_argument("--out", metavar="FILE", default="scree_plot.png", type=str,
                   help="The name of the output file [%(default)s]")


def safe_main():
    """A safe version of the main function (that catches ProgramError)."""
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Cancelled by user")
        sys.exit(0)
    except ProgramError as e:
        logger.error(e.message)
        parser.error(e.message)


# Calling the main, if necessary
if __name__ == "__main__":
    safe_main()
