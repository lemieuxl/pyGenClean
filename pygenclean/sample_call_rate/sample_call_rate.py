"""Removes samples with poor call rate."""


import logging
import argparse
from os import path

from ..utils.task import execute_external_command

from ..error import ProgramError

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "sample-call-rate"
DESCRIPTION = "Remove samples with poor call rate."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Remove samples with poor call rate.

    Args:
        args (argparse.Namespace): the arguments and options
        argv (list): the argument as list.

    These are the steps:

    1. Prints the options.
    2. Runs Plink with the ``mind`` option (:py:func:`runPlink`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Running Plink
    logger.info("Running Plink")
    run_plink(args)


def run_plink(options):
    """Run Plink with the ``mind`` option.

    Args:
        options (argparse.Namespace): the arguments and options.

    """
    # Executing the command
    execute_external_command(
        command=[
            "plink",
            "--noweb",
            "--bfile", options.bfile,
            "--mind", str(options.mind),
            "--make-bed",
            "--out", options.out,
        ]
    )


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the bed, bim and fam files
    for filename in [args.bfile + i for i in (".bed", ".bim", ".fam")]:
        if not path.isfile(filename):
            raise ProgramError(f"{filename}: no such file")

    # Check the mind option (between 0 and 1, inclusive)
    if (args.mind < 0) or (args.mind > 1):
        raise ProgramError(
            f"mind={args.mind}: must be between 0 and 1 (inclusive)",
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
        "--mind", type=float, metavar="FLOAT", default=0.1,
        help="The call rate threshold (remove samples with more than x "
             "percent missing genotypes). [Default: %(default)s]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="clean_mind",
        help="The prefix of the output files (wich will be a Plink binary "
             "file).  [default: %(default)s]",
    )
