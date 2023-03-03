"""Merges related samples according to IBS."""


import logging
import argparse
from os import path
from itertools import combinations
import random
from typing import List, Optional, Iterable, Dict, Tuple, Set

from ...error import ProgramError

from ...utils import timer

from ...version import pygenclean_version as __version__


SCRIPT_NAME = "merge-related-samples"
DESCRIPTION = "Merges related samples according to IBS."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Merges related samples according to IBS."""
    if args is None:
        args = parse_args(argv)
    check_args(args)

    merge_related_samples(args.ibs_related, args.out, args.no_status)


def merge_related_samples(filename: str, out_prefix: str,
                          no_status: bool) -> None:
    """Merge related samples.

    Args:
        filename (str): the name of the input file.
        out_prefix (str): the prefix of the output files.
        no_status (bool): is there a status column in the file?

    In the output file, there are a pair of samples per line. Hence, one can
    find related individuals by merging overlapping pairs.

    """
    # Generating the list of sample sets and status
    sample_sets = []
    sample_status = {}
    for row in file_parser(filename=filename, no_status=no_status):
        sample_sets.append({row["sample_1"], row["sample_2"]})

        if not no_status:
            sample_key = frozenset([row["sample_1"], row["sample_2"]])
            sample_status[sample_key] = row["status"]

    # Merging this list
    sample_sets = merge_sample_sets(sample_sets)

    chosen_samples = set()
    remaining_samples = set()

    # Printing the output file while selecting one samples for each group
    with open(out_prefix + ".merged_related_individuals", "w") as f:
        # Printing the header
        print("index", "FID1", "IID1", "FID2", "IID2", "status",
              sep="\t", file=f)

        for index, samples in enumerate(sample_sets):
            for two_samples in combinations(samples, 2):
                sample_key = frozenset(two_samples)
                status = sample_status.get(sample_key, "")
                print(index + 1, *two_samples[0], *two_samples[1], status,
                      sep="\t", file=f)

            # Choose a random sample from the group
            chosen = random.choice(list(samples))
            chosen_samples.add(chosen)
            remaining_samples |= samples - {chosen}

    # Printing the chosen samples
    with open(out_prefix + ".chosen_related_individuals", "w") as f:
        for sample_id in chosen_samples:
            print(*sample_id, sep="\t", file=f)

    # Printing the remaining samples
    with open(out_prefix + ".discarded_related_individuals", "w") as f:
        for sample_id in remaining_samples:
            print(*sample_id, sep="\t", file=f)


def merge_sample_sets(
    sample_sets: List[Set[Tuple[str, str]]],
) -> List[Set[Tuple[str, str]]]:
    """Merges the sample sets.

    From https://stackoverflow.com/a/9112588

    """
    was_merged = True
    while was_merged:
        was_merged = False
        results = []

        while sample_sets:
            common, rest = sample_sets[0], sample_sets[1:]
            sample_sets = []

            for x in rest:
                if x.isdisjoint(common):
                    sample_sets.append(x)

                else:
                    was_merged = True
                    common |= x

            results.append(common)

        sample_sets = results

    return sample_sets


def file_parser(filename: str,
                no_status: bool) -> Iterable[Dict[str, Tuple[str, str]]]:
    """Parses the file and returns the samples."""
    with open(filename) as f:
        header = None

        for line in f:
            row = line.rstrip().split("\t")

            if not header:
                header = {name: i for i, name in enumerate(row)}

                # Checking the column
                for col in ["FID1", "IID1", "FID2", "IID2"]:
                    if col not in header:
                        raise ProgramError(
                            f"{filename}: no column named '{col}'"
                        )

                # Checking the status
                if not no_status:
                    if "status" not in header:
                        raise ProgramError(
                            f"{filename}: no column named 'status'"
                        )

                continue

            to_yield = {
                "sample_1": (row[header["FID1"]], row[header["IID1"]]),
                "sample_2": (row[header["FID2"]], row[header["IID2"]]),
            }

            if not no_status:
                to_yield["status"] = row[header["status"]]

            yield to_yield


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not path.isfile(args.ibs_related):
        raise ProgramError(f"{args.ibs_related}: no such file")


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
    """Adds the arguments to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--ibs-related", type=str, metavar="FILE", required=True,
        help="The input file containing related individuals according to IBS "
             "value.",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--no-status", action="store_true",
        help="The input file doesn't have a 'status' column.",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="ibs_merged",
        help="The prefix of the output files. [default: %(default)s]",
    )
