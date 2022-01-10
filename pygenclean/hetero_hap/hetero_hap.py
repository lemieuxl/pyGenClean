"""Remove heterozygous haploid (males on chromosome X)."""


import logging
import argparse
from os import path

from ..utils.task import execute_external_command

from ..error import ProgramError

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "hetero-hap"
DESCRIPTION = "Remove heterozygous haploid (males on chromosome X)."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """The main function of this module.

    Args:
        args (argparse.Namespace): the arguments and options
        argv (list): the argument as a list.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Run plink
    logger.info("Running Plink to set heterozygous haploid as missing")
    run_plink(args)


def run_plink(options):
    """Sets heterozygous haploid markers to missing Plink.

    Args:
        options (argparse.Namespace): The arguments and options.

    """
    execute_external_command(
        command=[
            "plink",
            "--noweb",
            "--bfile", options.bfile,
            "--set-hh-missing",
            "--make-bed",
            "--out", options.out,
        ],
    )


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the tped and the tfam files
    for filename in [args.bfile + i for i in (".bed", ".bim", ".fam")]:
        if not path.isfile(filename):
            raise ProgramError(f"{filename}: no such file")


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
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively).",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="without_hh_genotypes",
        help="The prefix of the output files. [default: %(default)s]",
    )
