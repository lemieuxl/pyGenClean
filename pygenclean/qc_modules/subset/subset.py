"""Subset Plink data (samples or markers)."""


import argparse
import logging
from os import path
from typing import Dict, List, Optional, Tuple

from ...error import ProgramError
from ...utils import plink as plink_utils
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "subset"
DESCRIPTION = "Subset Plink datasets (samples or markers)."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Subset a Plink dataset (samples or markers).

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # The subset arguments
    samples, sample_subset_type = get_subset_arg(vars(args), "keep", "remove")
    markers, marker_subset_type = get_subset_arg(vars(args), "extract",
                                                 "exclude")

    # Executing
    logger.info("Subsetting data using Plink")
    plink_utils.subset(
        bfile=args.bfile,
        out=args.out,
        samples=samples,
        sample_subset_type=sample_subset_type,
        markers=markers,
        marker_subset_type=marker_subset_type,
        use_original_plink=args.plink_107,
    )


def get_subset_arg(
    kwargs: Dict[str, Optional[str]],
    first: str,
    second: str,
) -> Tuple[Optional[str], Optional[str]]:
    """Get the subset type."""
    for key in first, second:
        if kwargs[key]:
            return kwargs[key], key
    return None, None


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    Note:
        Only one operation for markers and one operation for samples can be
        executed at the same time.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    # Checking the subset files
    at_least_one = False
    for fn in (args.extract, args.exclude, args.remove, args.keep):
        if fn:
            at_least_one = True
            if not path.isfile(fn):
                raise ProgramError(f"{fn}: no such file")
    if not at_least_one:
        raise ProgramError("No operation was selected")


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Adds argument to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input files")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively).",
    )

    # The markers
    marker = group.add_mutually_exclusive_group()
    marker.add_argument(
        "--exclude", type=str, metavar="FILE",
        help="A file containing the list of markers to exclude from the "
             "dataset.",
    )
    marker.add_argument(
        "--extract", type=str, metavar="FILE",
        help="A file containing the list of markers to extract from the "
             "dataset.",
    )

    # The samples
    samples = group.add_mutually_exclusive_group()
    samples.add_argument(
        "--remove", type=str, metavar="FILE",
        help="A file containing the list of samples to remove from the "
             "dataset.",
    )
    samples.add_argument(
        "--keep", type=str, metavar="FILE",
        help="A file containing the list of samples to keep from the "
             "dataset.",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="subset",
        help="The prefix of the output files. [%(default)s]",
    )
