"""Compute heterozygosity rate and plots it."""


import argparse
import logging
from os import path
from typing import List, Optional, Tuple

import numpy as np
from geneparse.readers.plink import PlinkReader

from ...error import ProgramError
from ...utils import plink as plink_utils
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "heterozygosity-plot"
DESCRIPTION = "Computes heterozygosity rate and plots it."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None):
    """Compute heterozygosity rate and plots it.

    Args:
        args: The arguments and options
        argv: The argument as a list.

    These are the steps:

    1. Computes the heterozygosity rate.
    2. Save the heterozygosity rate for each samples in a file.
    3. Plot the heterozygosity rate (histogram or boxplot).

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
                        out_prefix: str):
    """Save the heterozygosity data to file.

    Args:
        heterozygosity: The heterozygosity data.
        samples: The list of samples.
        out_prefix: The prefix of the output files.

    """
    with open(out_prefix + ".het", "w") as f:
        for (fid, iid), het in zip(samples, heterozygosity):
            print(fid, iid, het, sep="\t", file=f)


def read_het_file(filename: str) -> np.ndarray:
    """Read the heterozygosity file.

    Args:
        filename: The name of the heterozygosity file.

    Returns:
        The heterozygosity rate for all samples.

    """
    data = []
    with open(filename) as f:
        for line in f:
            row = line.rstrip("\r\n").split()
            data.append(float(row[-1]))

    return np.array(data)


def compute_heterozygosity(prefix: str) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the heterozygosity ratio of all samples.

    Args:
        prefix: The prefix of the binary _Plink_ files.

    Returns:
        The heterozygosity rates for all the samples.

    """
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
                        options: argparse.Namespace):
    """Plot the heterozygosity rate distribution.

    Args:
        heterozygosity: The heterozygosity data.
        options: The options.

    Plot a histogram or a boxplot of the heterozygosity distribution.

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

    if options.xlim is not None:
        heterozygosity = heterozygosity[
            (heterozygosity >= options.xlim[0])
            & (heterozygosity <= options.xlim[1])
        ]

    # Plotting the histogram
    if options.boxplot:
        # Specific settings for the boxplot
        axe.yaxis.set_ticks_position("none")
        axe.spines["left"].set_visible(False)
        axe.yaxis.set_visible(False)
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
            weights=(
                np.zeros_like(heterozygosity) + 1.0 / heterozygosity.shape[0]
            ),
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
                                    lw=2, clip_on=False)
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


def check_args(args: argparse.Namespace):
    """Check the arguments and options.

    Args:
        args: The arguments and options.

    If there is a problem with an option, an exception is raised using the
    `ProgramError` class.

    """
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
    """Parse the command line options and arguments.

    Args:
        argv: An optional list of arguments.

    Returns:
        The parsed arguments and options.

    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser):
    """Add arguments and options to the parser.

    Args:
        parser: An argument parser to which arguments will be added.

    """
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
        help="The limit of the x axis (floats). The heterozygosity values"
             "will be restricted around those two values (inclusive).",
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
