"""Utility function related to Plink files."""


import re
from os import path

from . import decode_chrom, decode_sex
from .task import execute_external_command

from ..error import ProgramError


__all__ = ["check_files", "get_markers_on_chrom", "get_sample_sexes",
           "parse_bim", "parse_fam", "split_line", "extract_markers",
           "compare_bim"]


_SPACE_SPLITTER = re.compile(r"\s+")


def get_markers_on_chrom(bimfile, chromosomes):
    """Get the markers located on required chromosome(s) from the BIM file.

    Args:
        bimfile (str): the name of the BIM file.
        chromosomes (set): the list of chromosomes.

    Returns:
        set: the set of markers located on the required chromosomes.

    """
    return {
        marker.name for marker in parse_bim(bimfile)
        if marker.chrom in chromosomes
    }


def get_sample_sexes(famfile, only_iid=False):
    """Get sample's sex from the FAM file.

    Args:
        famfile (str): the name of the FAM file.

    Returns:
        dict: Sample to sex mapping.

    """
    sexes = {}
    for row in parse_fam(famfile):
        sample = row.iid if only_iid else (row.fid, row.iid)
        sexes[sample] = row.sex

    return sexes


def parse_bim(bimfile):
    """Parses a BIM file.

    Args:
        bimfile (str): the name of the BIM file.

    Returns:
        namedtuple: the row for each markers.

    """
    class BimRow():
        # pylint: disable=too-few-public-methods
        # pylint: disable=missing-docstring
        # pylint: disable=too-many-arguments
        # pylint: disable=invalid-name

        __slots__ = ["chrom", "name", "cm", "pos", "a1", "a2"]

        def __init__(self, chrom, name, cm, pos, a1, a2):
            self.chrom = decode_chrom(chrom)
            self.name = name
            self.cm = float(cm)
            self.pos = int(pos)
            self.a1 = a1
            self.a2 = a2

    with open(bimfile) as f:
        for line in f:
            yield BimRow(*line.rstrip("\r\n").split())


def parse_fam(famfile):
    """Parses a FAM file.

    Args:
        famfile (str): the name of the FAM file.

    Returns:
        object: the row for each sample.

    """
    class FamRow():
        # pylint: disable=too-few-public-methods
        # pylint: disable=too-many-arguments
        # pylint: disable=missing-docstring

        __slots__ = ["fid", "iid", "father", "mother", "sex", "status"]

        def __init__(self, fid, iid, father, mother, sex, status):
            self.fid = fid
            self.iid = iid
            self.father = father
            self.mother = mother
            self.sex = decode_sex(sex)
            self.status = status

        def __repr__(self):
            return f"<FamRow {self.fid}, {self.iid}>"

    with open(famfile) as f:
        for line in f:
            yield FamRow(*line.rstrip("\r\n").split())


def check_files(prefix):
    """Checks that all the files are present.

    Args:
        prefix (str): the prefix of the files.

    Returns:
        bool: ``True`` if all files are present, ``False`` otherwise.

    """
    for extension in ("bed", "bim", "fam"):
        if not path.isfile(f"{prefix}.{extension}"):
            return False
    return True


def split_line(line):
    """Split a Plink output line (spaces delimited).

    Args:
        line (str): the line to split.

    Returns:
        tuple: the split line.

    """
    return _SPACE_SPLITTER.split(line.strip())


def extract_markers(bfile, extract, out):
    """Extracts markers from a Plink file.

    Args:
        bfile (str): the prefix of the Plink file.
        extract (str): the name of the file containing the markers to extract.
        out (str): the prefix of the output file.

    """
    command = [
        "plink2",
        "--bfile", bfile,
        "--extract", extract,
        "--make-bed",
        "--out", out,
    ]
    execute_external_command(command)


def compare_bim(bim_a, bim_b):
    """Compares two BIM files.

    Args:
        bim_a (str): the first BIM file.
        bim_b (str): the second BIM file.

    Returns:
        tuple: A tuple of three sets. The first one contains the markers only in
        A. The second one contains the markers in both. The third one
        contains the markers only in B.

    """
    # Getting the set of markers in the the two BIM files
    markers_a = {row.name for row in parse_bim(bim_a)}
    markers_b = {row.name for row in parse_bim(bim_b)}

    # Markers in both BIM files
    in_both = markers_a & markers_b

    return (
        markers_a - in_both,
        in_both,
        markers_b - in_both,
    )
