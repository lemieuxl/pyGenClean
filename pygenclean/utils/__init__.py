"""Utility function which might be used throughout the pipeline."""


import gzip
import shlex
import functools
import subprocess

from ..error import ProgramError


__all__ = ["decode_chrom", "decode_sex", "is_gzip", "get_open_func",
           "execute_external_command", "split_extra_args"]


def decode_chrom(chrom):
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
        chrom = int(chrom)
    except ValueError:
        raise ProgramError(f"{chrom}: invalid chromosome")

    if chrom < 0 or chrom > 26:
        raise ProgramError(f"{chrom}: invalid chromosome")

    return chrom


def decode_sex(sex):
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


def is_gzip(filename):
    """Checks if a file is compressed using gzip.

    Args:
        filename (str): the file to check.

    Returns:
        bool: ``True`` if the file is compressed by gzip, false otherwise.

    """
    with open(filename, "rb") as f:
        first = f.read(2)
    return first == b'\x1f\x8b'


def get_open_func(filename):
    """Opens a file even if it's gzipped.

    Args:
        filename (str): the file to open

    Returns:
        file: the opened file

    """
    if is_gzip(filename):
        return functools.partial(gzip.open, filename, mode="rt")
    return functools.partial(open, filename)


def execute_external_command(*args):
    """Executes an external command.

    Args:
        args (list): the list containing the command and arguments/options.

    """

    try:
        proc = subprocess.Popen(
            args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        outs, errs = proc.communicate()

    except FileNotFoundError as exception:
        raise ProgramError(exception.strerror)

    if proc.returncode != 0:
        # Something went wrong
        command = " ".join(args)
        errs = errs.decode()
        raise ProgramError(f"Something went wrong:\n{errs}\n{command}")

    return outs


def split_extra_args(args):
    """Split extra arguments provided in a single string.

    Args:
        args (str): a string containing multiple arguments.

    Returns:
        list: the list of extra arguments (might be empty if args is None).

    """
    if args is None:
        return []
    return [shlex.quote(arg) for arg in shlex.split(args)]
