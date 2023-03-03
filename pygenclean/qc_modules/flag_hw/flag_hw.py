"""Flags markers failing Hardy Weinberg equilibrium."""


import logging
import argparse
from typing import List, Set, Optional

from ...utils.task import execute_external_command
from ...utils import plink as plink_utils
from ...utils import timer

from ...error import ProgramError

from ...version import pygenclean_version as __version__


SCRIPT_NAME = "flag-hw"
DESCRIPTION = "Flags markers failing Hardy Weinberg equilibrium."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Flags markers failing Hardy Weinberg equilibrium.

    Args:
        args (argparse.Namespace): the arguments and options
        argv (list): the argument as a list.

    These are the steps performed by this module:

    1. Prints the options of the module.
    2. Computes the number of markers in the input file
       (:py:func:`computeNumberOfMarkers`).
    3. If there are no markers, the module stops.
    4. Computes the Bonferroni therhold (:math:`0.05 / \\textrm{nbMarkers}`).
    5. Runs Plink to find failed markers with the Bonferroni threshold.
    6. Runs Plink to find failed markers with the default threshold.
    7. Compares the ``bim`` files for the Bonferroni threshold.
    8. Compares the ``bim`` files for the default threshold.
    9. Computes the "in between" marker list, which is the markers from the
       default threshold and the Bonferroni threshold.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Compute the number of markers
    logger.info("Counting the number of markers")
    nb_markers = get_nb_markers(args.bfile + ".bim")

    if nb_markers <= 0:
        logger.info("  - There are no markers: STOPPING NOW!")

    else:
        logger.info("  - There are %d markers", nb_markers)

        hw_threshold = str(0.05 / nb_markers)

        # Run the plink command for the Bonferonni threshold
        logger.info("Computing the HW equilibrium for %s", hw_threshold)
        compute_hw(args.bfile, hw_threshold,
                   args.out + ".threshold_" + hw_threshold, args.plink_107)

        # Run the plink command for the default threshold
        logger.info("Computing the HW equilibrium for %s", args.hwe)
        compute_hw(args.bfile, args.hwe, args.out + ".threshold_" + args.hwe,
                   args.plink_107)

        # Compare the BIM files
        logger.info("Creating the flagged SNP list for %s", hw_threshold)
        custom_snps = compare_bim(
            args.bfile + ".bim",
            args.out + ".threshold_" + hw_threshold + ".bim",
            args.out + ".snp_flag_threshold_" + hw_threshold,
        )

        logger.info("Creating the flagged SNP list for %s", args.hwe)
        hwe_snps = compare_bim(
            args.bfile + ".bim",
            args.out + ".threshold_" + args.hwe + ".bim",
            args.out + ".snp_flag_threshold_" + args.hwe,
        )

        logger.info("Creating the in between SNP list ([%s, %s[)",
                    args.hwe, hw_threshold)
        filename = (
            f"{args.out}.snp_flag_threshold_between_{args.hwe}-{hw_threshold}"
        )
        try:
            with open(filename, 'w') as f:
                differences = hwe_snps - custom_snps
                if len(differences) > 0:
                    print(*differences, sep="\n", file=f)

        except IOError as exception:
            raise ProgramError(f"{filename}: can't write file") from exception


def compare_bim(before: str, after: str, output_filename: str) -> Set[str]:
    """Compare two BIM files for differences.

    Args:
        before (str): the name of the file before modification.
        after (str): the name of the file after modification.
        output_filename: the name of the output file.

    Returns:
        set: the set difference before and after.

    The ``bim`` files contain the list of markers in a given dataset. The
    ``before`` file should have more markers than the ``after`` file. The
    ``after`` file should be a subset of the markers in the ``before`` file.

    """
    in_before, _, in_after = plink_utils.compare_bim(bim_a=before, bim_b=after)

    # There should be no extra markers in only 'after'
    if in_after:
        raise ProgramError(
            f"There should not be markers in {before} which are not "
            f"in {after}",
        )

    # We want to save the markers only in the first file
    with open(output_filename, "w") as f:
        print(*in_before, sep="\n", file=f)

    return in_before


def get_nb_markers(filename: str) -> int:
    """Count the number of markers (line) in a BIM file.

    Args:
        filename (str): the name of the BIM file.

    Returns:
        int: the number of markers (line) in the BIM file.

    """
    nb_line = 0
    with open(filename, "r") as f:
        nb_line = len(f.readlines())

    return nb_line


def compute_hw(prefix: str, threshold: str, out_prefix: str,
               use_original_plink: bool) -> None:
    """Compute the Hardy Weinberg test using Plink.

    Args:
        prefix (str): the prefix of all the files.
        threshold (str): the Hardy Weinberg threshold.
        out_prefix (str): the prefix of the output file.
        use_original_plink (bool): use plink 1.07 instead of plink 1.9.

    Uses Plink to exclude markers that failed the Hardy-Weinberg test at a
    specified significance threshold.

    """
    execute_external_command(
        command=[
            "plink" if use_original_plink else "plink1.9",
            "--noweb",
            "--bfile", prefix,
            "--hwe", threshold,
            "--make-bed",
            "--out", out_prefix,
        ]
    )


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :py:class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the required files
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    # Check the HW pvalue
    hw_value = args.hwe
    try:
        hw_value = float(hw_value)
    except ValueError as exception:
        raise ProgramError(f"{hw_value}: not a valid HW value") from exception

    if (hw_value < 0) or (hw_value > 1):
        raise ProgramError(f"{hw_value}: not a valid HW value")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}"
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
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
        "--hwe", type=str, metavar="FLOAT", default="1e-4",
        help="The Hardy-Weinberg equilibrium threshold. "
             "[default: %(default)s]",
    )
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="flag_hw",
        help="The prefix of the output files. [default: %(default)s]",
    )
