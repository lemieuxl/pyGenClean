"""Utility function related to Plink files."""


from shlex import quote

from . import decode_chrom, decode_sex
from ..utils import execute_external_command


__all__ = ["get_markers_on_chrom", "get_sample_sexes", "run_plink"]


def run_plink(*args):
    """Executes Plink with options.

    Args:
        args (list): the options for Plink.

    """
    command = ["plink", "--noweb"]
    for option in args:
        command.append(quote(option))
    execute_external_command(*command)


def get_markers_on_chrom(bimfile, chromosomes):
    """Get the markers located on required chromosome(s) from the BIM file.

    Args:
        bimfile (str): the name of the BIM file.
        chromosomes (set): the list of chromosomes.

    Returns:
        set: the set of markers located on the required chromosomes.

    """
    return {
        marker.name for marker in _parse_bim(bimfile)
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
    for row in _parse_fam(famfile):
        sample = row.iid if only_iid else (row.fid, row.iid)
        sexes[sample] = row.sex

    return sexes


def _parse_bim(bimfile):
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


def _parse_fam(famfile):
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

    with open(famfile) as f:
        for line in f:
            yield FamRow(*line.rstrip("\r\n").split())
