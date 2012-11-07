#!/usr/bin/env python2.7

import os
import sys
import gzip
import argparse

import numpy as npy

def main(argString=None):
    # Getting and checking the options
    args = parseArgs(argString)
    checkArgs(args)

    # Check the number of samples
    nb_samples = compute_nb_samples(args.tfile)

    # Compute the heterozygosity rate
    heterozygosity = compute_heterozygosity(args.tfile, nb_samples)

    # Plotting the heterozygosity rate distribution
    plot_heterozygosity(heterozygosity, args)


def compute_nb_samples(in_prefix):
    """Check the number of samples."""
    file_name = in_prefix + ".tfam"
    nb = None
    with open(file_name, 'rb') as input_file:
        nb = len(input_file.readlines())
    return nb


def is_heterozygous(genotype):
    """Tells if a genotype "A B" is heterozygous."""
    return not genotype[0] == genotype[-1]


def compute_heterozygosity(in_prefix, nb_samples):
    """Computes the heterozygosity ratio of samples (from tped)."""
    file_name = in_prefix + ".tped"

    # The function we want to use
    check_heterozygosity = npy.vectorize(is_heterozygous)

    # The autosomes
    autosomes = {str(i) for i in xrange(1, 23)}

    heterozygosity = npy.zeros(nb_samples, dtype=int)
    nb_markers = npy.zeros(nb_samples, dtype=int)
    with open(file_name, 'rb') as input_file:
        # There is no header
        for line in input_file:
            row = npy.array(line.rstrip("\r\n").split("\t"))

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

            if nb_markers[0] >= 10000:
                break

    return npy.true_divide(heterozygosity, nb_markers)


def plot_heterozygosity(heterozygosity, options):
    """Plots the heterozygosity rate distribution."""
    # importing important stuff
    import matplotlib as mpl
    if options.format != "X11":
        mpl.use("Agg")
    import matplotlib.pyplot as plt
    if options.format != "X11":
        plt.ioff()

    # Creating the figure and ax
    fig, ax = plt.subplots(1, 1, figsize=(13.5, 8))

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
    ax.set_ylabel("Proportion")

    # Plotting the histogram
    ax.hist(heterozygosity, bins=options.bins, color="#0099CC",
            histtype="stepfilled",  
            weights=npy.zeros_like(heterozygosity) + 1.0 / len(heterozygosity))

    # Plotting the mean, the median and the variance
    the_mean = npy.mean(heterozygosity)
    the_median = npy.median(heterozygosity)
    the_variance = npy.power(npy.std(heterozygosity), 2)
    mean_line = ax.axvline(the_mean, color="#CC0000", ls="--", lw=2,
                           clip_on=False)
    median_line = ax.axvline(the_median, color="#FF8800", ls="--", lw=2,
                             clip_on=False)
    variance_line = ax.axvline(the_variance, color="#FFFFFF", ls="--", lw=0,
                               clip_on=False)
    ax.legend([mean_line, median_line, variance_line],
              ["Mean ({:.4})".format(the_mean),
               "Median ({:.4})".format(the_median),
               "Variance ({:.4})".format(the_variance)],
              loc="best", prop={"size": 11})

    # The xlim
    if options.xlim is not None:
        ax.set_xlim(options.xlim)
    # The ylim
    if options.ymax is not None:
        ax.set_ylim(0.0, options.ymax)

    # Adding the mean of the heterozygosity

    # Saving the figure
    if options.format == "X11":
        plt.show()
    else:
        plt.savefig("{}.{}".format(options.out, options.format), dpi=300)

def checkArgs(args):
    """Checks the arguments and options.

    :param args: a :py:class:`Namespace` object containing the options of the
                 program.
    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed
    to the :class:`sys.stderr` and the program exists with code 1.

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


def parseArgs(argString=None): # pragma: no cover
    """Parses the command line options and arguments.

    :returns: A :py:class:`numpy.Namespace` object created by
              the :py:mod:`argparse` module. It contains the values of the
              different options.

    ===============  ======  ===================================================
        Options       Type                     Description
    ===============  ======  ===================================================
    ===============  ======  ===================================================

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
    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        return self.message

# The parser object
prog = "heterozygosity_plot"
desc = """Plot the distribution of the heterozygosity ratio."""
parser = argparse.ArgumentParser(description=desc, prog=prog)


# The INPUT files
group = parser.add_argument_group("Input File")
group.add_argument("--tfile", type=str, metavar="FILE", required=True,
                   help="The prefix of the transposed file")
# The options
group = parser.add_argument_group("Options")
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

args = parser.parse_args()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
