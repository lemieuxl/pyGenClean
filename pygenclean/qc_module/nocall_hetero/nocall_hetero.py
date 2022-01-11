"""Clean markers with no call or heterozygous only."""


import logging
import argparse
import shutil

import numpy as np

from geneparse.readers.plink import PlinkReader

from ...utils import plink as plink_utils

from ...error import ProgramError

from ...version import pygenclean_version as __version__


SCRIPT_NAME = "nocall-hetero"
DESCRIPTION = "Clean markers with no call or heterozygous only."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Clean markers with no call or heterozygous only.

    Args:
        args (argparse.Namespace): the arguments and options
        argv (list): the argument as a list.

    These are the steps:

    1. Prints the options.
    2. Reads the ``tfam`` and ``tped`` files and find all heterozygous and all
       failed markers (:py:func:`processTPEDandTFAM`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Process the files
    logger.info("Processing the binary Plink file")
    process_file(args.bfile, args.out)


def process_file(prefix, out_prefix):
    """Process the TPED and TFAM files.

    Args:
        prefix (str): the prefix of the BED/BIM/FAM files.
        out_prefix (str): the prefix of output files.

    Reads the input files and keeps in memory two sets containing the markers
    which are all failed or which contains only heterozygous genotypes.

    It creates two output files, ``prefix.allFailed`` and ``prefix.allHetero``,
    containing the markers that are all failed and are all heterozygous,
    respectively.

    .. note::
        All heterozygous markers located on the mitochondrial chromosome are
        not remove.

    """
    # The name of the bad SNPs
    all_failed = set()
    all_hetero = set()

    with PlinkReader(prefix) as bed:
        for data in bed.iter_genotypes():
            # Samples with no calls (i.e. NaN)
            no_calls = np.isnan(data.genotypes)

            # All genotyupes are 'no call'
            if no_calls.all():
                logger.debug("%s: all no call", data.variant.name)
                all_failed.add(data.variant.name)
                continue

            if data.variant.chrom != "MT":
                # Genotypes are either 'no call' or heterozygotes
                if ((data.genotypes == 1) | no_calls).all():
                    logger.debug("%s: all hetero", data.variant.name)
                    all_hetero.add(data.variant.name)
                    continue

    # Printing the SNPs with no calls
    with open(out_prefix + ".allFailed", 'w') as f:
        if len(all_failed) > 0:
            print(*all_failed, sep="\n", file=f)

    # Printing the SNPs with only hetero calls
    with open(out_prefix + ".allHetero", 'w') as f:
        if len(all_hetero) > 0:
            print(*all_hetero, sep="\n", file=f)

    # Printing the SNPs to exclude
    exclude = all_hetero | all_failed
    if len(exclude) > 0:
        logger.info("Excluding %d markers from original files", len(exclude))

        # Creating the exclude file
        filename = out_prefix + ".exclude"
        with open(filename, "w") as f:
            print(*exclude, sep="\n", file=f)

        # Excluding the markers
        plink_utils.subset_markers(prefix, filename, out_prefix, "exclude")

    else:
        logger.info("Copying original files (no markers to exclude)")
        for extension in (".bed", ".bim", ".fam"):
            shutil.copy(prefix + extension, out_prefix + extension)


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Checking the input files
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")


def parse_args(argv=None):
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser):
    """Add arguments and options to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the binary files by "
             "appending the prefix to .bed, .bim and .fam, respectively).",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="clean_noCall_hetero",
        help="The prefix of the output files. [default: %(default)s]",
    )
