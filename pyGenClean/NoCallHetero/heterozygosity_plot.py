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
import sys
import gzip
import logging
import argparse

import numpy as np

from .. import __version__


logger = logging.getLogger("heterozygosity_plot")


def main(argString=None):
    """The main function of the module.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. Checks the number of samples in the ``tfam`` file
       (:py:func:`compute_nb_samples`).
    3. Computes the heterozygosity rate (:py:func:`compute_heterozygosity`).
    4. Saves the heterozygosity data (in ``out.het``).
    5. Plots the heterozygosity rate (:py:func:`plot_heterozygosity`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    # Check the number of samples
    nb_samples = compute_nb_samples(args.tfile)

    # Compute the heterozygosity rate
    heterozygosity, samples = compute_heterozygosity(args.tfile, nb_samples)

    # Save heterozygosity data
    save_heterozygosity(heterozygosity, samples, args.out)

    # Plotting the heterozygosity rate distribution
    plot_heterozygosity(heterozygosity, args)


def save_heterozygosity(heterozygosity, samples, out_prefix):
    """Saves the heterozygosity data.

    :param heterozygosity: the heterozygosity data.
    :param samples: the list of samples.
    :param out_prefix: the prefix of the output files.

    :type heterozygosity: numpy.array
    :type samples: list of tuples of str
    :type out_prefix: str

    """
    # The output file
    o_file = None
    try:
        o_file = open(out_prefix + ".het", 'wb')
    except IOError:
        msg = "{}.het: can't write file".format(out_prefix)
        raise ProgramError(msg)

    # Saving the data
    for (fid, iid), het in zip(samples, heterozygosity):
        print >>o_file, "\t".join([fid, iid, str(het)])

    # Closing the file
    o_file.close()


def compute_nb_samples(in_prefix):
    """Check the number of samples.

    :param in_prefix: the prefix of the input file.

    :type in_prefix: str

    :returns: the number of sample in ``prefix.fam``.

    """
    file_name = in_prefix + ".tfam"
    nb = None
    with open(file_name, 'rb') as input_file:
        nb = len(input_file.readlines())
    return nb


def is_heterozygous(genotype):
    """Tells if a genotype "A B" is heterozygous.

    :param genotype: the genotype to test for heterozygosity.

    :type genotype: str

    :returns: ``True`` if the genotype is heterozygous, ``False`` otherwise.

    The genotype must contain two alleles, separated by a space. It then
    compares the first allele (``genotype[0]``) with the last one
    (``genotype[-1]``).

    .. testsetup::

        from pyGenClean.NoCallHetero.heterozygosity_plot import is_heterozygous

    .. doctest::

        >>> is_heterozygous("A A")
        False
        >>> is_heterozygous("G C")
        True
        >>> is_heterozygous("0 0") # No call is not heterozygous.
        False

    """
    return not genotype[0] == genotype[-1]


def compute_heterozygosity(in_prefix, nb_samples):
    """Computes the heterozygosity ratio of samples (from tped)."""
    tped_name = in_prefix + ".tped"
    tfam_name = in_prefix + ".tfam"

    # The function we want to use
    check_heterozygosity = np.vectorize(is_heterozygous)

    # The autosomes
    autosomes = {str(i) for i in xrange(1, 23)}

    # The tfam
    samples = None
    with open(tfam_name, 'rb') as input_file:
        samples = input_file.readlines()
        samples = [tuple(i.rstrip("\r\n").split("\t")[:2]) for i in samples]

    heterozygosity = np.zeros(nb_samples, dtype=int)
    nb_markers = np.zeros(nb_samples, dtype=int)
    with open(tped_name, 'rb') as input_file:
        # There is no header
        for line in input_file:
            row = np.array(line.rstrip("\r\n").split("\t"))

            chromosome = row[0]
            if chromosome not in autosomes:
                # This is not an autosome, so we skip
                continue

            # Getting the genotypes
            genotypes = row[4:]

            # Finding the heterozygous genotypes
            heterozygosity += check_heterozygosity(genotypes)

            # Adding to number of markers for each samples (excluding no calls)
            nb_markers += genotypes != "0 0"

    return np.true_divide(heterozygosity, nb_markers), samples


def plot_heterozygosity(heterozygosity, options):
    """Plots the heterozygosity rate distribution.

    :param heterozygosity: the heterozygosity data.
    :param options: the options.

    :type heterozygosity: numpy.array
    :type options: argparse.Namespace

    Plots a histogram or a boxplot of the heterozygosity distribution.

    """
    # importing important stuff
    import matplotlib as mpl
    if options.format != "X11" and mpl.get_backend() != "agg":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()

    # Creating the figure and ax
    figsize = (13.5, 8)
    if options.boxplot:
        figsize = (13.5, 4)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Setting subplot properties
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 9))
    ax.spines["left"].set_position(("outward", 9))

    # The title and labels
    ax.set_title(("Heterozygosity rate distribution "
                  "({:,} samples)".format(len(heterozygosity))), weight="bold")
    ax.set_xlabel("Heterozygosity rate")

    # Plotting the histogram
    if options.boxplot:
        # Specific settings for the boxplot
        ax.yaxis.set_ticks_position("none")
        ax.spines["left"].set_visible(False)
        ax.set_yticklabels([])
        fig.subplots_adjust(bottom=0.18)

        # The boxplot
        ax.boxplot(heterozygosity, notch=True, vert=False)
    else:
        # Specific settings for the histogram
        ax.set_ylabel("Proportion")

        # The histogram
        ax.hist(
            heterozygosity,
            bins=options.bins,
            color="#0099CC",
            histtype="stepfilled",
            weights=np.zeros_like(heterozygosity) + 1.0 / len(heterozygosity),
        )

        # Plotting the mean, the median and the variance
        the_mean = np.mean(heterozygosity)
        the_median = np.median(heterozygosity)
        the_variance = np.power(np.std(heterozygosity), 2)
        mean_line = ax.axvline(the_mean, color="#CC0000", ls="--", lw=2,
                               clip_on=False)
        median_line = ax.axvline(the_median, color="#FF8800", ls="--", lw=2,
                                 clip_on=False)
        variance_line = ax.axvline(the_variance, color="#FFFFFF", ls="--",
                                   lw=0, clip_on=False)
        ax.legend(
            [mean_line, median_line, variance_line],
            [
                "Mean ({:.4})".format(the_mean),
                "Median ({:.4})".format(the_median),
                "Variance ({:.4})".format(the_variance),
            ],
            loc="best",
            prop={"size": 11},
        )

        # The ylim
        if options.ymax is not None:
            ax.set_ylim(0.0, options.ymax)

    # The xlim
    if options.xlim is not None:
        ax.set_xlim(options.xlim)

    # Saving the figure
    if options.format == "X11":
        plt.show()
    else:
        file_name = "{}.{}".format(options.out, options.format)
        if options.boxplot:
            file_name = "{}_boxplot.{}".format(options.out, options.format)
        plt.savefig(file_name, dpi=300)


def checkArgs(args):
    """Checks the arguments and options.

    :param args: an object containing the options of the program.

    :type args: argparse.Namespace

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input file
    for file_name in [args.tfile + i for i in {".tfam", ".tped"}]:
        if not os.path.isfile(file_name):
            msg = "{}: no such file".format(file_name)
            raise ProgramError(msg)

    # Checking the xlimit
    if args.xlim is not None:
        if args.xlim[0] >= args.xlim[1]:
            msg = "invalid xlim"
            raise ProgramError(msg)

    # Checking the ymax
    if args.ymax is not None:
        if args.ymax <= 0:
            msg = "invalid ymax (must be above 0)"
            raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ============= ====== ======================================
       Options     Type              Description
    ============= ====== ======================================
    ``--tfile``   string The prefix of the transposed file.
    ``--boxplot`` bool   Draw a boxplot instead of a histogram.
    ``--format``  string The output file format.
    ``--bins``    int    The number of bins for the histogram.
    ``--xlim``    float  The limit of the x axis.
    ``--ymax``    float  "The maximal Y value.
    ``--out``     string The prefix of the output files.
    ============= ====== ======================================

    .. note::
        No option check is done here (except for the one automatically done by
        argparse). Those need to be done elsewhere (see :py:func:`checkArgs`).

    """
    args = None
    if argString is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argString)

    return args


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.

    :param msg: the message to print to the user before exiting.

    :type msg: str

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user

        :type msg: str

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


# The parser object
desc = "Plots the distribution of the heterozygosity ratio."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfile", type=str, metavar="FILE", required=True,
                   help="The prefix of the transposed file")
# The options
group = parser.add_argument_group("Options")
group.add_argument("--boxplot", action="store_true",
                   help=("Draw a boxplot instead of a histogram."))
group.add_argument("--format", type=str, metavar="FORMAT", default="png",
                   choices=["png", "ps", "pdf", "X11"],
                   help=("The output file format (png, ps, pdf, or X11 "
                         "formats are available). [default: %(default)s]"))
group.add_argument("--bins", type=int, metavar="INT", default=100,
                   help=("The number of bins for the histogram. [default: "
                         "%(default)d]"))
group.add_argument("--xlim", type=float, metavar="FLOAT", nargs=2,
                   help=("The limit of the x axis (floats)."))
group.add_argument("--ymax", type=float, metavar="FLOAT",
                   help=("The maximal Y value."))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="heterozygosity",
                   help=("The prefix of the output files. [default: "
                         "%(default)s]"))


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


if __name__ == "__main__":
    safe_main()
