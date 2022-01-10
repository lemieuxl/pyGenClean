"""Flag markers with MAF of 0 (monomorphic)."""


import logging
import argparse
from os import path

from ..utils import plink as plink_utils
from ..utils.task import execute_external_command

from ..error import ProgramError

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "flag-maf"
DESCRIPTION = "Flag markers with MAF of 0 (monomorphic)."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Flag monomorphic markers (i.e. MAF of 0).

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as list.

    These are the steps:

    1. Prints the options.
    2. Computes the frequencies using Plinl (:py:func:`computeFrequency`).
    3. Finds markers with MAF of 0, and saves them in a file
       (:py:func:`findSnpWithMaf0`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Compute frequency using plink
    logger.info("Computing the frequencies using Plink")
    compute_freq(args)

    # Read the freqency file
    logger.info("Flagging SNPs with MAF = 0")
    find_maf_0(args.out + ".frq", args.out)


def find_maf_0(freq_filename, prefix):
    """Finds SNPs with MAF of 0 and put them in a file.

    Args:
        freq_filename (str): the name of the frequency file.
        prefix (str): the prefix of all the files.

    Reads a frequency file from Plink, and find markers with a minor allele
    frequency of zero.

    """
    if not path.isfile(freq_filename):
        raise ProgramError(f"{freq_filename}: no such file")

    maf_0_set = set()
    na_set = set()

    with open(freq_filename, "r") as f:
        header = None
        for line in f:
            row = plink_utils.split_line(line)

            if header is None:
                header = {name: i for i, name in enumerate(row)}
                for col_name in ("SNP", "MAF"):
                    if col_name not in header:
                        raise ProgramError(
                            f"{freq_filename}: no column named {col_name}",
                        )

                continue

            # We have data
            name = row[header["SNP"]]
            maf = row[header["MAF"]]

            if maf == "0":
                # We want to flag this SNP
                maf_0_set.add(name)

            elif maf == "NA":
                # We want to flag this SNP, because the MAF is NA
                na_set.add(name)

    # Creating the output files
    if len(maf_0_set) == 0:
        logger.info("  - There are no markers with MAF 0")
    else:
        logger.info("  - There are %s markers with MAF 0", len(maf_0_set))

    with open(prefix + ".list", "w") as output_file:
        for marker_name in maf_0_set:
            print(marker_name, file=output_file)

    if len(na_set) > 0:
        logger.info("  - There are %s markers with NA MAF", len(na_set))
        with open(prefix + ".na_list", "w") as output_file:
            for marker_name in na_set:
                print(marker_name,  file=output_file)


def compute_freq(options):
    """Compute the frequency of the SNPs.

    Args:
        options (argparse.Namespace): the arguments and options.

    """
    # Executing the command
    execute_external_command(
        command=[
            "plink",
            "--noweb",
            "--bfile", options.bfile,
            "--freq",
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


def parse_args(argv=None):
    """Parses the arguments and function."""
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
             "respectively.",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="flag_maf_0",
        help="The prefix of the output files. [default: %(default)s]",
    )