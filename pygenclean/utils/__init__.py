"""Utility function which might be used throughout the pipeline."""


import datetime
import functools
import gzip
import shlex
import time
from typing import Any, Callable, Iterable, List, Optional, Set

from ..error import ProgramError


__all__ = ["decode_chrom", "decode_sex", "is_gzip", "get_open_func",
           "split_extra_args", "timer", "flip_alleles", "compatible_alleles",
           "check_allele_status"]


def decode_chrom(chrom: str) -> int:
    """Decode chromosomes from string to integer.

    Args:
        chrom (str): the chromosome to encode.

    Returns:
        int: the decoded chromosome.

    It changes ``X``, ``Y``, ``XY`` and ``MT`` to ``23``, ``24``, ``25`` and
    ``26``, respectively. It changes everything else as :py:class:`int`.

    An exception is raised if :py:class:`ValueError` is raised, or if the
    decoded chromosome is below 0 or above 26.

    Note:
        We allow for chromosome 0, as Illumina uses these to "remove" a marker.

    """
    chrom = chrom.upper()

    if chrom == "X":
        return 23

    if chrom == "Y":
        return 24

    if chrom == "XY":
        return 25

    if chrom == "MT":
        return 26

    try:
        num_chrom = int(chrom)
    except ValueError as exception:
        raise ProgramError(f"{chrom}: invalid chromosome") from exception

    if num_chrom < 0 or num_chrom > 26:
        raise ProgramError(f"{chrom}: invalid chromosome")

    return num_chrom


def decode_sex(sex: str) -> str:
    """Decodes the sex of a sample.

    Args:
        sex (str): the sex to decode.

    Returns:
        str: the decoded sex.

    Notes:
        It changes ``1`` and ``2`` to ``Male`` and ``Female`` respectively. It
        encodes everything else to ``Unknown``.

    """
    if sex == "1":
        return "Male"

    if sex == "2":
        return "Female"

    return "Unknown"


def is_gzip(filename: str) -> bool:
    """Checks if a file is compressed using gzip.

    Args:
        filename (str): the file to check.

    Returns:
        bool: ``True`` if the file is compressed by gzip, false otherwise.

    """
    with open(filename, "rb") as f:
        first = f.read(2)
    return first == b'\x1f\x8b'


def get_open_func(filename: str) -> functools.partial:
    """Opens a file even if it's gzipped.

    Args:
        filename (str): the file to open

    Returns:
        file: the opened file

    """
    if is_gzip(filename):
        return functools.partial(gzip.open, filename, mode="rt")
    return functools.partial(open, filename)


def split_extra_args(args: str) -> List[str]:
    """Split extra arguments provided in a single string.

    Args:
        args (str): a string containing multiple arguments.

    Returns:
        list: the list of extra arguments (might be empty if args is None).

    """
    if args is None:
        return []
    return list(map(shlex.quote, shlex.split(args)))


def timer(decorated_logger):
    """A simple timer decorator which logs the time."""
    def inner(function: Callable) -> Callable:
        """A simple timer decorator."""
        @functools.wraps(function)
        def wrapper_timer(*args: list, **kwargs: dict) -> Any:
            """The wrapper for the decorator."""
            # Timing
            start_time = time.perf_counter()
            value = function(*args, **kwargs)
            end_time = time.perf_counter()

            # Logging
            run_time = datetime.timedelta(seconds=end_time - start_time)
            decorated_logger.info("Done in %s", run_time)

            # Returning the original value
            return value

        return wrapper_timer

    return inner


_COMPLEMENT_ALLELE = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "0": "0",
}


def flip_alleles(alleles: Iterable[str]) -> Set[str]:
    """Flip alleles in a set."""
    try:
        return {_COMPLEMENT_ALLELE[allele] for allele in alleles}

    except KeyError as error:
        raise ProgramError(f"{alleles}: unknown alleles") from error


def compatible_alleles(alleles_1: Set[str], alleles_2: Set[str]) -> bool:
    """Check if two sets of alelels are compatible."""
    as_is = len((alleles_1 | alleles_2) - {"0"}) <= 2
    as_comp = len((alleles_1 | flip_alleles(alleles_2)) - {"0"}) <= 2
    return as_is or as_comp


def check_allele_status(alleles_1: Set[str],
                        alleles_2: Set[str]) -> Optional[str]:
    """Check the status between two sets of alleles."""
    # Removing the missing allele (i.e. "0")
    alleles_1 -= {"0"}
    alleles_2 -= {"0"}

    # The final status
    status = None

    # Different alleles
    if alleles_1 != alleles_2:
        # Same number of alleles, trying to flip one set
        if len(alleles_1) == len(alleles_2):
            alleles_2 = flip_alleles(alleles_2)
            if alleles_1 == alleles_2:
                if len(alleles_1) == 1:
                    status = "homo_flip"
                else:
                    status = "flip"
            else:
                status = "problem"

        # Different number of alleles
        else:
            # One allele in common
            if len(alleles_1 & alleles_2) == 1:
                status = "homo_hetero"
            else:
                # Trying to flip one set
                alleles_2 = flip_alleles(alleles_2)
                if len(alleles_1 & alleles_2) == 1:
                    status = "homo_hetero_flip"
                else:
                    status = "problem"

    return status
