"""Flag markers failing Hardy Weinberg equilibrium."""


import argparse
import logging
from typing import Dict, List, Optional, Set

from ...error import ProgramError
from ...report.summaries import FlagHwSummary
from ...utils import plink as plink_utils
from ...utils import timer
from ...utils.task import execute_external_command
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "flag-hw"
DESCRIPTION = "Flags markers failing Hardy Weinberg equilibrium."
DEFAULT_OUT = "flag_hw"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> Dict[str, Optional[str]]:
    """Flag markers failing Hardy Weinberg equilibrium.

    Args:
        args: the arguments and options
        argv: the argument as a list.

    Returns:
        A dictionary containing summary information about the run.

    These are the steps performed by this module:

    1. Compute the number of markers in the input file.
    2. Compute the Bonferroni therhold ($0.05 / \\textrm{nb markers}$).
    3. Run _Plink_ to find failed markers with the Bonferroni threshold.
    4. Run _Plink_ to find failed markers with the default threshold.
    5. Compare the `bim` files for the Bonferroni threshold.
    6. Compare the `bim` files for the default threshold.
    7. Compute the "in between" marker list, which is the markers from the
       default threshold and the Bonferroni threshold.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    logger.info("%s", DESCRIPTION)

    # Compute the number of markers
    logger.info("Counting the number of markers")
    nb_markers = get_nb_markers(args.bfile + ".bim")

    # The summary
    summary = FlagHwSummary(args)

    if nb_markers <= 0:
        logger.info("  - There are no markers: STOPPING NOW!")
        return {
            "summary": summary,
            "usable_files": {
                "bfile": args.bfile,
            },
        }

    logger.info("  - There are %d markers", nb_markers)

    hw_threshold = str(0.05 / nb_markers)

    # Run the plink command for the Bonferroni threshold
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

    return {
        "summary": summary,
        "usable_files": {
            "bfile": args.bfile,
            "flagged": args.out + ".snp_flag_threshold_" + args.hwe,
            "flagged_bonferroni": args.out + ".snp_flag_threshold_"
                                  + hw_threshold,
        },
    }


def compare_bim(before: str, after: str, output_filename: str) -> Set[str]:
    """Compare two BIM files for differences.

    Args:
        before: The name of the file before modification.
        after: The name of the file after modification.
        output_filename: The name of the output file.

    Returns:
        The set difference before and after.

    The `bim` files contain the list of markers in a given dataset. The
    `before` file should have more markers than the `after` file. The
    `after` file should be a subset of the markers in the `before` file.

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
        filename: The name of the BIM file.

    Returns:
        The number of markers (line) in the BIM file.

    """
    nb_line = 0
    with open(filename, "r") as f:
        nb_line = len(f.readlines())

    return nb_line


def compute_hw(prefix: str, threshold: str, out_prefix: str,
               use_original_plink: bool):
    """Compute the Hardy Weinberg test using Plink.

    Args:
        prefix: The prefix of all the files.
        threshold: The Hardy Weinberg threshold.
        out_prefix: The prefix of the output file.
        use_original_plink: Use the original Plink (version 1.07).

    Use Plink to exclude markers that failed the Hardy-Weinberg test at a
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


def check_args(args: argparse.Namespace):
    """Check the arguments and options.

    Args:
        args: The arguments and options.

    If there is a problem with an option, an exception is raised using the
    `ProgramError` class.

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
    """Parse the command line options and arguments.

    Args:
        argv: An optional list of arguments.

    Returns:
        The parsed arguments and options.

    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}"
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
        "--out", type=str, metavar="FILE", default=DEFAULT_OUT,
        help="The prefix of the output files. [default: %(default)s]",
    )
