"""Utility function related to Plink files."""


import re
from os import path
from typing import Tuple, Set, List

from . import decode_chrom, decode_sex
from .task import execute_external_command


__all__ = ["check_files", "get_markers_on_chrom", "get_sample_sexes",
           "parse_bim", "parse_fam", "split_line", "subset_markers",
           "compare_bim"]


_SPACE_SPLITTER = re.compile(r"\s+")


def get_markers_on_chrom(bimfile: str, chromosomes: Set[int]) -> Set[str]:
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


def get_sample_sexes(famfile: str, only_iid: bool = False) -> dict:
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


class BimRow():
    # pylint: disable=too-few-public-methods
    # pylint: disable=missing-docstring
    # pylint: disable=too-many-arguments
    # pylint: disable=invalid-name

    __slots__ = ["chrom", "name", "cm", "pos", "a1", "a2"]

    def __init__(self, chrom: str, name: str, cm: str, pos: str, a1: str,
                 a2: str):
        self.chrom = decode_chrom(chrom)
        self.name = name
        self.cm = float(cm)
        self.pos = int(pos)
        self.a1 = a1
        self.a2 = a2

    def __repr__(self) -> str:
        return f"<BimRow {self.name}>"


def parse_bim(bimfile: str) -> BimRow:
    """Parses a BIM file.

    Args:
        bimfile (str): the name of the BIM file.

    Returns:
        namedtuple: the row for each markers.

    """
    with open(bimfile) as f:
        for line in f:
            yield BimRow(*line.rstrip("\r\n").split())


class FamRow():
    # pylint: disable=too-few-public-methods
    # pylint: disable=too-many-arguments
    # pylint: disable=missing-docstring

    __slots__ = ["fid", "iid", "father", "mother", "sex", "status"]

    def __init__(self, fid: str, iid: str, father: str, mother: str, sex: str,
                 status: str):
        self.fid = fid
        self.iid = iid
        self.father = father
        self.mother = mother
        self.sex = decode_sex(sex)
        self.status = status

    def __repr__(self) -> str:
        return f"<FamRow {self.fid}, {self.iid}>"


def parse_fam(famfile: str) -> FamRow:
    """Parses a FAM file.

    Args:
        famfile (str): the name of the FAM file.

    Returns:
        object: the row for each sample.

    """
    with open(famfile) as f:
        for line in f:
            yield FamRow(*line.rstrip("\r\n").split())


def check_files(prefix: str) -> bool:
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


def split_line(line: str) -> List[str]:
    """Split a Plink output line (spaces delimited).

    Args:
        line (str): the line to split.

    Returns:
        list: the split line.

    """
    return _SPACE_SPLITTER.split(line.strip())


def subset_markers(bfile: str, markers: str, out: str, subset_type: str,
                   use_original_plink: bool = False) -> None:
    """Extracts markers from a Plink file.

    Args:
        bfile (str): the prefix of the Plink file.
        markers (str): the name of the file containing the markers to subset.
        out (str): the prefix of the output file.
        subset_type (str): either `extract` or `exclude`.

    """
    if subset_type not in {"exclude", "extract"}:
        raise ValueError(f"{subset_type}: invalid subset")

    command = [
        "plink" if use_original_plink else "plink1.9",
        "--noweb",
        "--bfile", bfile,
        f"--{subset_type}", markers,
        "--make-bed",
        "--out", out,
    ]
    execute_external_command(command)


def compare_bim(bim_a: str, bim_b: str) -> Tuple[Set[str], Set[str], Set[str]]:
    """Compares two BIM files.

    Args:
        bim_a (str): the first BIM file.
        bim_b (str): the second BIM file.

    Returns:
        tuple: A tuple of three sets. The first one contains the markers only
               in A. The second one contains the markers in both. The third one
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
