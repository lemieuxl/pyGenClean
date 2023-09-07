"""Clean markers with no call or heterozygous genotypes only."""


import argparse
import logging
import shutil
from typing import Dict, List, Optional

import numpy as np
from geneparse.readers.plink import PlinkReader

from ...error import ProgramError
from ...report.summaries import NoCallHeteroSummary
from ...utils import plink as plink_utils
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "nocall-hetero"
DESCRIPTION = "Clean markers with no call or heterozygous only."
DEFAULT_OUT = "clean_noCall_hetero"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> Dict[str, str]:
    """Clean markers with no call or heterozygous genotypes only.

    Args:
        args: The arguments and options
        argv: The argument as a list.

    Returns:
        A dictionary containing summary information about the run.

    These are the steps:

    1. Read the binary _Plink_ files and find all heterozygous and all failed
       markers.
    2. Use _Plink_ to exclude the markers found at the first step.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    logger.info("%s", DESCRIPTION)

    # Process the files
    process_file(args.bfile, args.out, args.plink_107)

    return {
        "summary": NoCallHeteroSummary(args),
        "usable_files": {
            "bfile": args.out,
        },
    }


def process_file(prefix: str, out_prefix: str, use_original_plink: bool):
    """Process the binary _Plink_ files.

    Args:
        prefix: The prefix of the BED/BIM/FAM files.
        out_prefix: The prefix of the output files.

    Reads the input files and keeps in memory two sets containing the markers
    which are all failed or which contains only heterozygous genotypes.

    It creates two output files, `{prefix}.all_failed` and
    `{prefix}.all_hetero`, containing the markers that are all failed and are
    all heterozygous, respectively.

    note:
        All heterozygous markers located on the mitochondrial chromosome are
        not removed.

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
    with open(out_prefix + ".all_failed", 'w') as f:
        if len(all_failed) > 0:
            print(*all_failed, sep="\n", file=f)

    # Printing the SNPs with only hetero calls
    with open(out_prefix + ".all_hetero", 'w') as f:
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
        plink_utils.subset(
            bfile=prefix,
            out=out_prefix,
            markers=filename,
            marker_subset_type="exclude",
            use_original_plink=use_original_plink,
        )

    else:
        logger.info("Copying original files (no markers to exclude)")
        for extension in (".bed", ".bim", ".fam"):
            shutil.copy(prefix + extension, out_prefix + extension)


def check_args(args: argparse.Namespace):
    """Check the arguments and options.

    Args:
        args: The arguments and options.

    If there is a problem with an option, an exception is raised using the
    `ProgramError` class.

    """
    # Checking the input files
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")


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
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the binary files by "
             "appending the prefix to .bed, .bim and .fam, respectively).",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default=DEFAULT_OUT,
        help="The prefix of the output files. [default: %(default)s]",
    )
