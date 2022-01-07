"""Removes markers with poor call rate."""


import logging
import argparse
from os import path

from ..utils.task import execute_external_command
from ..utils import plink as plink_utils

from ..error import ProgramError

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "marker-call-rate"
DESCRIPTION = "Remove markers with poor call rate."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Removes markers with poor call rate.

    Args:
        args (arparse.Namespace): the arguments and options
        argv (list): the argument as a list.

    These are the steps:

    1. Prints the options.
    2. Runs Plink with the ``geno`` option (:py:func:`runPlink`).
    3. Compares the two ``bim`` files (before and after the Plink ``geno``
       analysis) (:py:func:`compareBIM`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Run plink
    logger.info("Running Plink")
    run_plink(args)

    # Comparing the bim
    logger.info("Comparing BIM files")
    compare_bim(args)


def compare_bim(args):
    """Compare two BIM files.

    Args:
        args (argparse.Namespace): the arguments and options.

    Compares two BIM files (before and after a change).

    """
    # Comparing the two BIM files
    in_before, _, in_after = plink_utils.compare_bim(
        bim_a=args.bfile + ".bim",
        bim_b=args.out + ".bim",
    )

    # There should be no extra markers in only 'after'
    if in_after:
        raise ProgramError(
            f"There should not be markers in {args.bfile}.bim which are not "
            f"in {args.out}.bim",
        )

    # We want to save the the markers only in the first file
    with open(args.out + ".removed_snps", "w") as f:
        print(*in_before, sep="\n", file=f)


def run_plink(options):
    """Runs Plink with the ``geno`` option.

    Args:
        options (argparse.Namespace): the arguments and options.

    """
    # Executing the command
    execute_external_command(
        command=[
            "plink",
            "--noweb",
            "--bfile", options.bfile,
            "--geno", str(options.geno),
            "--make-bed",
            "--out",
            options.out,
        ]
    )


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the tped and the tfam files
    for filename in [args.bfile + i for i in (".bed", ".bim", ".fam")]:
        if not path.isfile(filename):
            raise ProgramError(f"{filename}: no such file")

    # Check the mind option (between 0 and 1, inclusive)
    if (args.geno < 0) or (args.geno > 1):
        raise ProgramError(
            f"geno={args.geno}: must be between 0 and 1 (inclusive)",
        )


def parse_args(argv=None):
    """Parses the arguments and options."""
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
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively).",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--geno", type=float, metavar="FLOAT", default=0.02,
        help="The missingness threshold (remove SNPs with more than x percent "
             "missing genotypes). [Default: %(default)s]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="clean_geno",
        help="The prefix of the output files. [default: %(default)s]",
    )
