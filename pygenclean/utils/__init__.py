"""Utility function which might be used throughout the pipeline."""


import gzip
import time
import datetime
import shlex
import functools
from typing import Iterable, List, Any, Callable, Set

from ..error import ProgramError


__all__ = ["decode_chrom", "decode_sex", "is_gzip", "get_open_func",
           "split_extra_args", "timer", "flip_alleles"]


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
    return [shlex.quote(arg) for arg in shlex.split(args)]


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
