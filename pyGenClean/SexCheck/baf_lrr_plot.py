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


logger = logging.getLogger("baf_llr_plot")


def main(argString=None):
    """The main function of this module.

    :param argString: the options.

    :type argString: list

    These are the steps:

    1. Prints the options.
    2. Reads the problematic samples (:py:func:`read_problematic_samples`).
    3. Finds and checks the raw files for each of the problematic samples
       (:py:func:`check_file_names`).
    4. Plots the BAF and LRR (:py:func:`plot_baf_lrr`).

    """
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    # Reading the problematic samples
    samples = read_problematic_samples(args.problematic_samples)

    # Checking for file names
    file_names = check_file_names(samples, args.raw_dir, args)

    # Plotting the BAF and LRR for chromosome X and Y
    plot_baf_lrr(file_names, args)


def check_file_names(samples, raw_dir, options):
    """Check if all files are present.

    :param samples: a list of tuples with the family ID as first element (str)
                    and sample ID as last element (str).
    :param raw_dir: the directory containing the raw files.
    :param options: the options.

    :type samples: list of tuples
    :type raw_dir: str
    :type options: argparse.Namespace

    :returns: a dict containing samples as key (a tuple with the family ID as
              first element and sample ID as last element) and the name of the
              raw file as element.

    """
    file_names = {}
    for sample in samples:
        the_sample = None
        try:
            the_sample = sample[1]
        except IndexError:
            msg = ("problematic samples file should include both family and "
                   "individual IDs")
            raise ProgramError(msg)
        if options.use_full_ids:
            the_sample = options.full_ids_delimiter.join(sample)

        file_name = os.path.join(raw_dir, "{}.txt".format(the_sample))
        if not os.path.isfile(file_name):
            file_name += ".gz"
            if not os.path.isfile(file_name):
                msg = "can't find file for sample {}".format(the_sample)
                raise ProgramError(msg)
        file_names[the_sample] = file_name

    return file_names


def read_problematic_samples(file_name):
    """Reads a file with sample IDs.

    :param file_name: the name of the file containing problematic samples after
                      sex check.

    :type file_name: str

    :returns: a set of problematic samples (tuple containing the family ID as
              first element and the sample ID as last element).

    Reads a file containing problematic samples after sex check. The file is
    provided by the module :py:mod:`pyGenClean.SexCheck.sex_check`. This file
    contains two columns, the first one being the family ID and the second one,
    the sample ID.

    """
    problematic_samples = set()
    open_func = open
    if file_name.endswith(".gz"):
        open_func = gzip.open
    with open_func(file_name, 'rb') as input_file:
        problematic_samples = set([
            tuple(i.rstrip("\r\n").split("\t")) for i in input_file.readlines()
        ])

    return problematic_samples


def encode_chromosome(chromosome):
    """Encodes chromosomes.

    :param chromosome: the chromosome to encode.

    :type chromosome: str

    :returns: the encoded chromosome.

    Encodes the sexual chromosomes, from ``23`` and ``24`` to ``X`` and ``Y``,
    respectively.

    .. note::
        Only the sexual chromosomes are encoded.

    .. testsetup::

        from pyGenClean.SexCheck.baf_lrr_plot import encode_chromosome

    .. doctest::

        >>> encode_chromosome("23")
        'X'
        >>> encode_chromosome("24")
        'Y'
        >>> encode_chromosome("This is not a chromosome")
        'This is not a chromosome'

    """
    if chromosome == "23":
        return "X"
    if chromosome == "24":
        return "Y"
    return chromosome


def plot_baf_lrr(file_names, options):
    """Plot BAF and LRR for a list of files.

    :param file_names: contains the name of the input file for each sample.
    :param options: the options.

    :type file_names: dict
    :type options: argparse.Namespace

    Plots the BAF (B Allele Frequency) and LRR (Log R Ratio) of each samples.
    Only the sexual chromosome are shown.

    """
    # importing important stuff
    import matplotlib as mpl
    if options.format != "X11" and mpl.get_backend() != "agg":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()

    # For each of the sample/files
    for sample, file_name in file_names.iteritems():
        data = []

        # Reading the file
        open_func = open
        if file_name.endswith(".gz"):
            open_func = gzip.open
        with open_func(file_name, 'rb') as input_file:
            header_index = dict([
                (col_name, i)
                for i, col_name in
                enumerate(input_file.readline().rstrip("\r\n").split("\t"))
            ])
            for col_name in {"Chr", "Position", "B Allele Freq",
                             "Log R Ratio"}:
                if col_name not in header_index:
                    msg = "{}: no column named {}".format(file_name, col_name)
                    raise ProgramError(msg)

            # Reading the dat
            for line in input_file:
                row = line.rstrip("\r\n").split("\t")

                # We only need X and Y chromosomes
                chromosome = encode_chromosome(row[header_index["Chr"]])
                if chromosome not in {"X", "Y"}:
                    continue

                # The position
                position = row[header_index["Position"]]
                try:
                    position = int(position)
                except ValueError:
                    msg = "{}: impossible position {}".format(file_name,
                                                              position)
                    raise ProgramError(msg)

                # The BAF
                baf = row[header_index["B Allele Freq"]]
                try:
                    baf = float(baf)
                except ValueError:
                    msg = "{}: impossible baf {}".format(file_name, baf)
                    raise ProgramError(msg)

                # The LRR
                lrr = row[header_index["Log R Ratio"]]
                try:
                    lrr = float(lrr)
                except ValueError:
                    msg = "{}: impossible lrr {}".format(file_name, lrr)
                    raise ProgramError(msg)

                # Saving the data
                data.append((chromosome, position, lrr, baf))

        # Creating the numpy array
        data = np.array(data, dtype=[("chr", "a1"), ("pos", int),
                                     ("lrr", float), ("baf", float)])

        # Creating the figure and axes
        fig, axes = plt.subplots(2, 2, figsize=(20, 8))
        plt.subplots_adjust(left=0.05, right=0.97, wspace=0.15, hspace=0.3)
        fig.suptitle(sample, fontsize=16, weight="bold")

        # Setting subplot properties
        for ax in axes.flatten():
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_position(("outward", 9))
            ax.spines["left"].set_position(("outward", 9))

        # Separating the axes
        x_lrr_ax, x_baf_ax, y_lrr_ax, y_baf_ax = axes.flatten(order='F')

        # Printing the X chromosome
        curr_chr = data["chr"] == "X"
        x_lrr_ax.plot(data["pos"][curr_chr]/1000000.0, data["lrr"][curr_chr],
                      "o", ms=1, mec="#0099CC",
                      mfc="#0099CC")[0].set_clip_on(False)
        x_baf_ax.plot(data["pos"][curr_chr]/1000000.0, data["baf"][curr_chr],
                      "o", ms=1, mec="#669900",
                      mfc="#669900")[0].set_clip_on(False)
        x_lrr_ax.axhline(y=0, color="#000000", ls="--", lw=1.2)
        x_baf_ax.axhline(y=0.5, color="#000000", ls="--", lw=1.2)
        x_lrr_ax.set_ylabel("LRR", weight="bold")
        x_baf_ax.set_ylabel("BAF", weight="bold")
        x_baf_ax.set_xlabel("Position (Mb)", weight="bold")
        x_lrr_ax.set_title("Chromosome X", weight="bold")

        # Printing the X chromosome
        curr_chr = data["chr"] == "Y"
        y_lrr_ax.plot(data["pos"][curr_chr]/1000000.0, data["lrr"][curr_chr],
                      "o", ms=1, mec="#0099CC",
                      mfc="#0099CC")[0].set_clip_on(False)
        y_baf_ax.plot(data["pos"][curr_chr]/1000000.0, data["baf"][curr_chr],
                      "o", ms=1, mec="#669900",
                      mfc="#669900")[0].set_clip_on(False)
        y_lrr_ax.axhline(y=0, color="#000000", ls="--", lw=1.2)
        y_baf_ax.axhline(y=0.5, color="#000000", ls="--", lw=1.2)
        y_lrr_ax.set_ylabel("LRR", weight="bold")
        y_baf_ax.set_ylabel("BAF", weight="bold")
        y_baf_ax.set_xlabel("Position (Mb)", weight="bold")
        y_lrr_ax.set_title("Chromosome Y", weight="bold")

        # Saving the figure
        if options.format == "X11":
            plt.show()
        else:
            plt.savefig(
                "{}_{}_lrr_baf.{}".format(options.out, sample, options.format),
                dpi=options.dpi,
            )

        # Closing the figure
        plt.close(fig)


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
    if not os.path.isfile(args.problematic_samples):
        msg = "{}: no such file".format(args.problematic_samples)
        raise ProgramError(msg)

    # Checking the raw directory
    if not os.path.isdir(args.raw_dir):
        msg = "{}: no such directory".format(args.raw_dir)
        raise ProgramError(msg)

    # Checking the DPI value
    if args.dpi < 10:
        msg = "{}: DPI too low".format(args.dpi)
        raise ProgramError(msg)

    return True


def parseArgs(argString=None):  # pragma: no cover
    """Parses the command line options and arguments.

    :param argString: the options.

    :type argString: list

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the
              different options.

    ========================= ====== ==========================================
              Options          Type                     Description
    ========================= ====== ==========================================
    ``--problematic-samples`` string The list of sample with sex problems to
                                     plot
    ``--use-full-ids``        bool   Use full sample IDs (famID and indID).
    ``--full-ids-delimiter``  string The delimiter between famID and indID.
    ``--raw-dir``             string Directory containing information about
                                     every samples (BAF and LRR).
    ``--format``              string The output file format (png, ps, pdf, or
                                     X11).
    ``--out``                 string The prefix of the output files.
    ========================= ====== ==========================================

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
desc = "Plots the BAF and LRR of problematic samples."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("-v", "--version", action="version",
                    version="pyGenClean version {}".format(__version__))

# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--problematic-samples", type=str, metavar="FILE",
                   required=True,
                   help=("A file containing the list of samples with sex "
                         "problems (family and individual ID required, "
                         "separated by a single tabulation). Uses "
                         "only the individual ID by default, unless "
                         "--use-full-ids is used."))
group.add_argument("--use-full-ids", action="store_true",
                   help=("Use this options to use full sample IDs (famID and "
                         "indID). Otherwise, only the indID will be use."))
group.add_argument("--full-ids-delimiter", type=str, metavar="CHAR",
                   default="_", help=("The delimiter between famID and indID "
                                      "for the raw file names. [default: "
                                      "%(default)s]"))
group.add_argument("--raw-dir", type=str, metavar="DIR", required=True,
                   help=("Directory containing information about every "
                         "samples (BAF and LRR)."))
# The options
group = parser.add_argument_group("Options")
group.add_argument("--format", type=str, metavar="FORMAT", default="png",
                   choices=["png", "ps", "pdf", "X11"],
                   help=("The output file format (png, ps, pdf, or X11 "
                         "formats are available). [default: %(default)s]"))
group.add_argument("--dpi", type=int, metavar="DPI", default=300,
                   help=("The pixel density of the figure(s) (DPI). "
                         "[default: %(default)d]"))
# The OUTPUT files
group = parser.add_argument_group("Output File")
group.add_argument("--out", type=str, metavar="FILE",
                   default="sexcheck",
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
