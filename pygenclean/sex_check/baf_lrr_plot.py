"""Plots the BAF and LRR of samples with sex mismatch."""


import logging
import argparse
from os import path
from glob import glob

import pandas as pd

import matplotlib.pyplot as plt

from ..utils import decode_chrom
from ..utils import illumina

from ..error import ProgramError

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "baf-lrr-plot"
DESCRIPTION = "Plots the BAF and LRR of samples with sex mismatch."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Plots the BAF and LRR of samples with sex mismatch.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1. Reads the samples with sex mismatch
       (:py:func:`read_problematic_samples`).
    2. Finds and checks the raw files for each of the samples
       (:py:func:`check_file_names`).
    3. Plots the BAF and LRR (:py:func:`plot_baf_lrr`).

    Note:
        Since the intensity file only contain a single sample ID (instead of
        the traditional ``FID`` and ``IID``, only the ``IID`` is kept.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reading the sample with sex mismatches
    logger.info("Getting samples with sex mismatch")
    samples = read_samples(args.problematic_samples)

    # Searching for the intensities
    logger.info("Searching for intensity files")
    filenames = get_intensity_filenames(samples, args.intensity_dir, args)

    # Plotting the BAF and LRR values
    logger.info("Plotting LRR and BAF values")
    for sample, filename in filenames.items():
        plot_baf_lrr(sample, filename, args)


def plot_baf_lrr(sample, filename, args):
    """Plots BAF and LRR values for multiple samples.

    Args:
        sample (str): the sample ID.
        filename (str): the name of the input file.
        args (argparse.Namespace): the arguments and options.

    """
    # pylint: disable=no-member
    # The columns
    cols = {
        "chrom": args.intensity_chrom_col,
        "pos": args.intensity_pos_col,
        "baf": args.intensity_baf_col,
        "lrr": args.intensity_lrr_col,
    }

    # Checking the first hundred lines to check if [Data] tag is present
    start = illumina.get_data_position(filename)

    # Reading the file (there might be a [Data] flag)
    df = pd.read_csv(
        filename, sep=args.intensity_delimiter, skiprows=start,
        usecols=[cols["chrom"], cols["pos"], cols["baf"], cols["lrr"]],
        dtype={cols["chrom"]: str}
    )

    # Decoding the chromosomes
    df["__chrom"] = df.loc[:, cols["chrom"]].map(decode_chrom)

    # Keeping only sexual chromosomes
    df = df.loc[df.loc[:, "__chrom"].isin({23, 24}), :]

    # Plotting
    figure, axes = plt.subplots(2, 2, figsize=(20, 8))
    plt.subplots_adjust(wspace=0.15, hspace=0.3)
    figure.suptitle(sample, fontsize=16, weight="bold")

    # Setting subplot properties
    for axe in axes.flatten():
        axe.xaxis.set_ticks_position("bottom")
        axe.yaxis.set_ticks_position("left")
        axe.spines["top"].set_visible(False)
        axe.spines["right"].set_visible(False)
        axe.spines["bottom"].set_position(("outward", 9))
        axe.spines["left"].set_position(("outward", 9))

    # The chromosomes
    chrx = df.loc[:, "__chrom"] == 23
    chry = df.loc[:, "__chrom"] == 24

    # LRR
    axe = axes[0, 0]
    axe.scatter(
        df.loc[chrx, cols["pos"]] / 1e6, df.loc[chrx, cols["lrr"]], s=1,
        c="#0099CC",
    )
    axe.set_title("Chromosome X", weight="bold")
    axe.set_ylabel("LRR", weight="bold")

    # BAF
    axe = axes[1, 0]
    axe.scatter(
        df.loc[chrx, cols["pos"]] / 1e6, df.loc[chrx, cols["baf"]], s=1,
        c="#669900",
    )
    axe.set_xlabel("Position (Mb)", weight="bold")
    axe.set_ylabel("BAF", weight="bold")

    # LRR
    axe = axes[0, 1]
    axe.scatter(
        df.loc[chry, cols["pos"]] / 1e6, df.loc[chry, cols["lrr"]], s=1,
        c="#0099CC",
    )
    axe.set_title("Chromosome Y", weight="bold")
    axe.set_ylabel("LRR", weight="bold")

    # BAF
    axe = axes[1, 1]
    axe.scatter(
        df.loc[chry, cols["pos"]] / 1e6, df.loc[chry, cols["baf"]], s=1,
        c="#669900",
    )
    axe.set_xlabel("Position (Mb)", weight="bold")
    axe.set_ylabel("BAF", weight="bold")

    # Saving the figure
    plt.savefig(f"{args.out}_{sample}_lrr_baf.{args.format}", dpi=args.dpi,
                bbox_inches="tight")
    plt.close(figure)


def read_samples(filename):
    """Reads the samples from a file.

    Args:
        filename (str): the name of the file.

    Returns:
        set: a set of problematic samples (tuple, first element is FID and
        second is IID.

    Reads a file containing problematic samples after sex-check. The file is
    provided by the module :py:mod:`pygenclean.sex_check`. This file contains
    two columns (white-space delimited), the first one being FID and the second
    one, IID.

    """
    with open(filename) as f:
        samples = {tuple(line.split()) for line in f.read().splitlines()}
    return samples


def get_intensity_filenames(samples, intensity_dir, args):
    """Searches for intensity file for each sample.

    Args:
        samples (str): the samples (FID, IID).
        intensity_dir (str): the directory containing the intensity files.
        args (argparse.Namespace): the arguments and options.

    Returns:
        dict: Samples as key (a tuple with the family ID as first element and
        sample ID as last element) and the name of the intensity file.

    """
    filenames = {}
    for fid, iid in samples:
        sample = iid
        if args.use_full_ids:
            sample = fid + args.full_ids_delimiter + iid

        filename = glob(path.join(
            intensity_dir, f"{sample}.*",
        ))

        if len(filename) == 0:
            raise ProgramError(f"{sample}: no intensity file found")

        if len(filename) != 1:
            raise ProgramError(f"{sample}: multiple intensity files found")

        filenames[sample] = filename.pop()

    return filenames


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not path.isfile(args.problematic_samples):
        raise ProgramError(f"{args.problematic_samples}: no such file")

    # Checking the raw directory
    if not path.isdir(args.intensity_dir):
        raise ProgramError(f"{args.intensity_dir}: no such directory")

    # Checking the DPI value
    if args.dpi < 100:
        raise ProgramError(f"{args.dpi}: DPI too low")


def parse_args(argv=None):
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments
    add_args(parser)

    if argv is None:
        return parser.parse_args()
    return parser.parse_args(argv)


def add_args(parser):
    """Adds argument to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input files")
    group.add_argument(
        "--problematic-samples", type=str, metavar="FILE", required=True,
        help="A file containing the list of samples with sex problems (family "
             "and individual ID required, separated by a white space). "
             "Uses only the individual ID by default, unless "
             "``--use-full-ids`` is used.",
    )

    # The intensity file options
    group = parser.add_argument_group("Intensity file options")
    group.add_argument(
        "--use-full-ids", action="store_true",
        help="Use this options to use full sample IDs (famID and indID). "
             "Otherwise, only the indID will be use.",
    )
    group.add_argument(
        "--full-ids-delimiter", type=str, metavar="CHAR", default="_",
        help="The delimiter between famID and indID for the intensity file "
             "names. [%(default)s]",
    )
    group.add_argument(
        "--intensity-dir", type=str, metavar="DIR", required=True,
        help="Directory containing information about every samples (BAF and "
             "LRR).",
    )
    group.add_argument(
        "--intensity-delimiter", type=str, metavar="SEP", default="\t",
        help="The field separator for the intensity file. [tabulation]",
    )
    group.add_argument(
        "--intensity-chrom-col", type=str, metavar="COL", default="Chr",
        help="The name of the column containing chromosomes in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-pos-col", type=str, metavar="COL", default="MapInfo",
        help="The name of the column containing positions in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-baf-col", type=str, metavar="COL",
        default="B Allele Freq",
        help="The name of the column containing BAF in the intensity "
             "file. [%(default)s]",
    )
    group.add_argument(
        "--intensity-lrr-col", type=str, metavar="COL",
        default="Log R Ratio",
        help="The name of the column containing LRR in the intensity "
             "file. [%(default)s]",
    )

    # The options
    group = parser.add_argument_group("Plot options")
    group.add_argument(
        "--format", type=str, metavar="FORMAT", default="png",
        choices=["png", "ps", "pdf"],
        help="The output file format (png, ps or pdf formats are available). "
             "[%(default)s]",
    )
    group.add_argument(
        "--dpi", type=int, metavar="DPI", default=300,
        help="The pixel density of the figure(s) (DPI). [%(default)d]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="sexcheck",
        help="The prefix of the output files. [%(default)s]",
    )
