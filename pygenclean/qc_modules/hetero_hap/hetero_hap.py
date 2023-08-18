"""Remove heterozygous haploid genotypes (males on chromosome X)."""


import argparse
import logging
from typing import Dict, List, Optional

from ...error import ProgramError
from ...report.summaries import HeteroHapSummary
from ...utils import plink as plink_utils
from ...utils import timer
from ...utils.task import execute_external_command
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "hetero-hap"
DESCRIPTION = "Remove heterozygous haploid (males on chromosome X)."
DEFAULT_OUT = "without_hh_genotypes"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> Dict[str, str]:
    """Remove heterozygous haploid genotypes (males on chromosome X).

    Args:
        args: the arguments and options
        argv: the argument as a list.

    Returns:
        A dictionary containing summary information about the run.

    These are the steps:

    1. Use Plink to remove heterozygous haploid genotypes.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Logging
    logger.info("Running Plink to set heterozygous haploid as missing")
    logger.info("  --bfile '%s'", args.bfile)
    logger.info("  --out '%s'", args.out)

    run_plink(args)

    # Returns a dictionary of usable files (for next step, if any)
    return {
        "summary": HeteroHapSummary(args),
        "usable_files": {
            "bfile": args.out,
        },
    }


def run_plink(options: argparse.Namespace):
    """Set heterozygous haploid markers to missing Plink.

    Args:
        options: The arguments and options.

    """
    execute_external_command(
        command=[
            "plink" if options.plink_107 else "plink1.9",
            "--noweb",
            "--bfile", options.bfile,
            "--set-hh-missing",
            "--make-bed",
            "--out", options.out,
        ],
    )


def check_args(args: argparse.Namespace):
    """Check the arguments and options.

    Args:
        args: The arguments and options.

    If there is a problem with an option, an exception is raised using the
    `ProgramError` class.

    """
    # Check if we have the tped and the tfam files
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
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively).",
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
