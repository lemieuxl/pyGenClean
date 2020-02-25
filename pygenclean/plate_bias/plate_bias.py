"""Finds plate bias (if any)."""


import logging
import argparse
from os import path

from ..error import ProgramError

from ..utils import plink as plink_utils
from ..utils.task import execute_external_commands

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "plate-bias"
DESCRIPTION = "Check for plate bias."


logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Plots the BAF and LRR of samples with sex mismatch.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1. Runs a plate bias analysis using Plink.
    2. Extracts the list of significant markers after plate bias analysis.
    3. Computes the frequency of all significant markers after plate bias
       analysis.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reading the plates
    sample_plates = get_plates(args.plates)

    # Creating the plate files
    create_plate_files(sample_plates, args.bfile + ".fam", args.out)

    # Executing Plink on each plate file
    execute_plate_bias(args.bfile, set(sample_plates.values()), args.p_filter,
                       args.out, args.nb_threads)


def execute_plate_bias(bfile, plates, p_filter, out, threads):
    """Execute plate bias on each plates.

    Args:

    """
    commands = [
        ["plink", "--noweb", "--bfile", bfile, "--fisher",
         "--pheno", f"{out}.{plate}.pheno", "--pfilter", str(p_filter),
         "--out", f"{out}.{plate}"]
        for plate in plates
    ]
    logger.info("Executing %d analysis", len(commands))
    execute_external_commands(commands, threads=threads)


def create_plate_files(sample_plates, fam, prefix):
    """Creates the files for the test (one file per plate.

    Args:
        sample_plates (dict): the plate for each sample.
        fam (str): the fam file.
        prefix (str): the output prefix.

    """
    logger.info("Generating plate bias files")
    samples = list(plink_utils.parse_fam(fam))

    for plate in set(sample_plates.values()):
        with open(f"{prefix}.{plate}.pheno", "w") as f:
            for sample in samples:
                sample_plate = sample_plates[(sample.fid, sample.iid)]
                status = 2 if sample_plate == plate else 1

                print(sample.fid, sample.iid, status, file=f)


def get_plates(filename):
    """Gets the plate for each samples.

    Args:
        filename (str): the file containing the plates.

    Returns:
        dict: the (FID, IID) assigned to each plate.

    """
    logger.info("Reading plate information from '%s'", filename)
    plates = {}
    with open(filename) as f:
        for line in f:
            fid, iid, plate = line.rstrip().split()
            plates[(fid, iid)] = plate
    return plates


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: missing plink files")

    if not path.isfile(args.plates):
        raise ProgramError(f"{args.plates}: no such file")


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
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively.",
    )
    group.add_argument(
        "--plates", type=str, metavar="FILE", required=True,
        help="The file containing the plate organization of each samples. "
             "Must contains three column (with no header): famID, indID and "
             "plateName.",
    )
    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--p-filter", type=float, metavar="FLOAT", default=1e-7,
        help="The significance threshold used for the plate effect. "
             "[%(default).1e]",
    )
    group.add_argument(
        "--nb-threads", type=int, metavar="N", default=1,
        help="The number of threads for this analysis. [%(default)d]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="plate_bias",
        help="The prefix of the output files. [%(default)s]",
    )
