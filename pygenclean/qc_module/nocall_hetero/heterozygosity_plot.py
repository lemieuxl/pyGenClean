"""Computes heterozygosity rate and plots it."""


import logging
import argparse
from os import path
from typing import Optional, List, Tuple

import numpy as np

from geneparse.readers.plink import PlinkReader

from ...utils import plink as plink_utils

from ...error import ProgramError
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "heterozygosity-plot"
DESCRIPTION = "Computes heterozygosity rate and plots it."


logger = logging.getLogger(__name__)


def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Computes heterozygosity rate and plots it.

    Args:
        args (argparse.Namespace): the arguments and options
        argv (list): the argument as a list.

    These are the steps:

    1. Prints the options.
    2. Checks the number of samples in the ``tfam`` file
       (:py:func:`compute_nb_samples`).
    3. Computes the heterozygosity rate (:py:func:`compute_heterozygosity`).
    4. Saves the heterozygosity data (in ``out.het``).
    5. Plots the heterozygosity rate (:py:func:`plot_heterozygosity`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    heterozygosity = None
    if args.het_file is not None:
        heterozygosity = read_het_file(args.het_file)

    else:
        # Compute the heterozygosity rate and save it to file
        heterozygosity, samples = compute_heterozygosity(args.bfile)
        save_heterozygosity(heterozygosity, samples, args.out)

    # Plotting the heterozygosity rate distribution
    plot_heterozygosity(heterozygosity, args)


def save_heterozygosity(heterozygosity: np.ndarray, samples: np.ndarray,
                        out_prefix: str) -> None:
    """Saves the heterozygosity data.

    Args:
        heterozygosity (numpy.array): the heterozygosity data.
        samples (list): the list of samples.
        out_prefix (str): the prefix of the output files.

    """
    with open(out_prefix + ".het", "w") as f:
        for (fid, iid), het in zip(samples, heterozygosity):
            print(fid, iid, het, sep="\t", file=f)


def read_het_file(filename: str) -> np.ndarray:
    """Reads the heterozygosity file.

    Args:
        filename (str): the name of the heterozygosity file.

    Returns:
        numpy.ndarray: the heterozygosity of all samples.

    """
    data = []
    with open(filename) as f:
        for line in f:
            row = line.rstrip("\r\n").split()
            data.append(float(row[-1]))

    return np.array(data)


def compute_heterozygosity(prefix: str) -> Tuple[np.ndarray, np.ndarray]:
    """Computes the heterozygosity ratio of samples."""
    logger.info("Computing heterozygosity from binary files")
    # The autosomes
    autosomes = {str(i) for i in range(1, 23)}

    # Parsing the binary files
    with PlinkReader(prefix) as bed:
        fam = bed.fam

        # Extracting the samples
        samples = fam.reset_index().loc[:, ["fid", "iid"]].to_numpy()

        # The number of heterygous and total markers (zeros)
        nb_hetero = np.zeros(fam.shape[0], dtype=int)
        nb_markers = np.zeros(fam.shape[0], dtype=int)

        for data in bed.iter_genotypes():
            # Keeping only the autosomes
            if data.variant.chrom not in autosomes:
                continue

            # Counting the number of heterozygous and all markers
            nb_hetero[data.genotypes == 1] += 1
            nb_markers[~np.isnan(data.genotypes)] += 1

    return np.true_divide(nb_hetero, nb_markers), samples


def plot_heterozygosity(heterozygosity: np.ndarray,
                        options: argparse.Namespace) -> None:
    """Plots the heterozygosity rate distribution.

    Args:
        heterozygosity (numpy.ndarray): the heterozygosity data.
        options (argparse.Namespace): the options.

    Plots a histogram or a boxplot of the heterozygosity distribution.

    """
    # pylint: disable=import-outside-toplevel

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
    figure, axe = plt.subplots(1, 1, figsize=figsize)

    # Setting subplot properties
    axe.xaxis.set_ticks_position("bottom")
    axe.yaxis.set_ticks_position("left")
    axe.spines["top"].set_visible(False)
    axe.spines["right"].set_visible(False)
    axe.spines["bottom"].set_position(("outward", 9))
    axe.spines["left"].set_position(("outward", 9))

    # The title and labels
    axe.set_title(
        f"Heterozygosity rate distribution ({heterozygosity.shape[0]:,} "
        "samples)",
        weight="bold",
    )
    axe.set_xlabel("Heterozygosity rate")

    # Plotting the histogram
    if options.boxplot:
        # Specific settings for the boxplot
        axe.yaxis.set_ticks_position("none")
        axe.spines["left"].set_visible(False)
        axe.set_yticklabels([])
        figure.subplots_adjust(bottom=0.18)

        # The boxplot
        axe.boxplot(heterozygosity, notch=True, vert=False)

    else:
        # Specific settings for the histogram
        axe.set_ylabel("Proportion")

        # The histogram
        axe.hist(
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
        mean_line = axe.axvline(the_mean, color="#CC0000", ls="--", lw=2,
                                clip_on=False)
        median_line = axe.axvline(the_median, color="#FF8800", ls="--", lw=2,
                                  clip_on=False)
        variance_line = axe.axvline(the_variance, color="#FFFFFF", ls="--",
                                    lw=0, clip_on=False)
        axe.legend(
            [mean_line, median_line, variance_line],
            [
                f"Mean ({the_mean:.4})",
                f"Median ({the_median:.4})",
                f"Variance ({the_variance:.4})",
            ],
            loc="best",
            prop={"size": 11},
        )

        # The ylim
        if options.ymax is not None:
            axe.set_ylim(0.0, options.ymax)

    # The xlim
    if options.xlim is not None:
        axe.set_xlim(options.xlim)

    # Saving the figure
    if options.format == "X11":
        plt.show()

    else:
        file_name = f"{options.out}.{options.format}"
        if options.boxplot:
            file_name = f"{options.out}_boxplot.{options.format}"
        plt.savefig(file_name, dpi=300)


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options."""
    if args.het_file is None and args.bfile is None:
        raise ProgramError("Specify either '--het-file' or '--bfile'")

    if args.het_file is not None:
        if not path.isfile(args.het_file):
            raise ProgramError(f"{args.het_file}: no such file")

    # Checking the input file
    if args.bfile is not None:
        if not plink_utils.check_files(args.bfile):
            raise ProgramError(f"{args.bfike}: no such binary files")

    # Checking the xlimit
    if args.xlim is not None:
        if args.xlim[0] >= args.xlim[1]:
            raise ProgramError("invalid xlim")

    # Checking the ymax
    if args.ymax is not None:
        if args.ymax <= 0:
            raise ProgramError("invalid ymax (must be above 0)")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments and options to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--bfile", type=str, metavar="FILE",
        help="The prefix of the transposed file",
    )
    group.add_argument(
        "--het-file", type=str, metavar="FILE",
        help="The heterozygosity file created by this QC module.",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--boxplot", action="store_true",
        help="Draw a boxplot instead of a histogram.",
    )
    group.add_argument(
        "--format", type=str, metavar="FORMAT", default="png",
        choices=["png", "ps", "pdf", "X11"],
        help="The output file format (png, ps, pdf, or X11 formats are "
             "available). [default: %(default)s]",
    )
    group.add_argument(
        "--bins", type=int, metavar="INT", default=100,
        help="The number of bins for the histogram. [default: %(default)d]",
    )
    group.add_argument(
        "--xlim", type=float, metavar="FLOAT", nargs=2,
        help="The limit of the x axis (floats).",
    )
    group.add_argument(
        "--ymax", type=float, metavar="FLOAT",
        help="The maximal Y value.",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="heterozygosity",
        help="The prefix of the output files. [default: %(default)s]",
    )
