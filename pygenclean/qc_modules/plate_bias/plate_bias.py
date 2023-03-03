"""Finds plate bias (if any)."""


import re
import logging
import argparse
from os import path
from glob import glob
from typing import Optional, List, Dict, Tuple, Set

from ...error import ProgramError

from ...utils import task
from ...utils import plink as plink_utils
from ...utils import timer

from ...version import pygenclean_version as __version__


SCRIPT_NAME = "plate-bias"
DESCRIPTION = "Check for plate bias."


logger = logging.getLogger(__name__)


class SignificantMarker():
    # pylint: disable=too-few-public-methods
    # pylint: disable=missing-docstring
    # pylint: disable=too-many-arguments
    # pylint: disable=invalid-name

    __slots__ = ["chrom", "pos", "name", "p_value", "odds", "plate"]

    def __init__(self, chrom: str, pos: str, name: str, p_value: str,
                 odds: str, plate: str):
        self.chrom = chrom
        self.pos = pos
        self.name = name
        self.p_value = p_value
        self.odds = odds
        self.plate = plate


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Finds plate biases (if any).

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1. Runs a plate bias analysis using Plink.
    2. Extracts the list of significant markers after plate bias analysis.
    3. Computes the frequency of all significant markers after plate bias
       analysis.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Reading the plates
    sample_plates = get_plates(args.plates)

    # Creating the plate files
    create_plate_files(sample_plates, args.bfile + ".fam", args.out)

    # Executing Plink on each plate file
    execute_plate_bias(args.bfile, set(sample_plates.values()), args.p_filter,
                       args.out, args.nb_tasks, args.plink_107)

    # Extracting significant markers
    assoc_results = extract_significant_markers(args.out)

    # Computing the frequencies of the significant markers
    maf = compute_frequency(args.bfile, args.out, args.plink_107)

    write_summary(assoc_results, maf, args.out)


def write_summary(results: List[SignificantMarker], maf: Dict[str, str],
                  out: str) -> None:
    """Writes the results.

    Args:
        results (list): the plate bias results.
        maf (dict): the MAF of each markers.
        out (str): the output rpefix.

    """
    with open(out + ".significant_markers.summary.tsv", "w") as f:
        print("chrom", "pos", "name", "maf", "p", "odds", "plate",
              sep="\t", file=f)

        for marker in results:
            print(marker.chrom, marker.pos, marker.name, maf[marker.name],
                  marker.p_value, marker.odds, marker.plate, sep="\t", file=f)


def compute_frequency(bfile: str, out: str,
                      use_original_plink: bool = False) -> Dict[str, str]:
    """Compute the frequency of specific markers.

    Args:
        bfile (str): the prefix of the input Plink file.
        out (str): the prefix of the output files.

    Returns:
        dict: the MAF for each marker.

    """
    logger.info("Computing the frequencies of significant markers")

    plink_utils.compute_freq(bfile=bfile, out=out + ".significant_markers",
                             use_original_plink=use_original_plink)

    # Reading the MAF
    frequencies = {}
    with open(out + ".significant_markers.frq", "r") as f:
        header = None
        for line in f:
            row = plink_utils.split_line(line)

            if header is None:
                header = {name: i for i, name in enumerate(row)}
                continue

            marker = row[header["SNP"]]
            maf = row[header["MAF"]]

            frequencies[marker] = maf

    return frequencies


def extract_significant_markers(prefix: str) -> List[SignificantMarker]:
    """Extract significant markers from the association files.

    Args:
        prefix (str): the prefix of the association files.

    Returns:
        list: the list of significant markers.

    """
    data = []

    assoc_files = glob(prefix + ".*.assoc.fisher")

    for filename in assoc_files:
        # Finding plate name
        plate_name = re.search(r"plate_bias\.(\S+)\.assoc\.fisher$", filename)
        if plate_name is None:
            raise ProgramError(f"{filename}: impossible to find plate name")
        plate_name = plate_name.group(1)

        # Reading the file
        with open(filename) as f:
            header = None
            for line in f:
                row = plink_utils.split_line(line)

                if header is None:
                    header = {name: i for i, name in enumerate(row)}
                    continue

                data.append(SignificantMarker(
                    chrom=row[header["CHR"]],
                    pos=row[header["BP"]],
                    name=row[header["SNP"]],
                    p_value=row[header["P"]],
                    odds=row[header["OR"]],
                    plate=plate_name,
                ))

    # Creating the output file
    with open(prefix + ".significant_markers.txt", "w") as f:
        print(*{marker.name for marker in data}, sep="\n", file=f)

    return data


def execute_plate_bias(bfile: str, plates: Set[str], p_filter: float,
                       out: str, threads: int,
                       use_original_plink: bool = False) -> None:
    """Execute plate bias on each plates.

    Args:
        bfile (str): the file prefix.
        plates (set): the set of plates.
        p_filter (float): the p-value filter.
        out (str): the output prefix.
        threads (int): the number of threads.

    """
    commands = [
        [
            "plink" if use_original_plink else "plink1.9",
            "--noweb",
            "--bfile", bfile,
            "--fisher",
            "--pheno", f"{out}.{plate}.pheno",
            "--pfilter", str(p_filter),
            "--out", f"{out}.{plate}",
        ] for plate in plates
    ]
    logger.info("Executing %d analysis", len(commands))
    task.execute_external_commands(commands, threads=threads)


def create_plate_files(sample_plates: Dict[Tuple[str, str], str],
                       fam: str, prefix: str) -> None:
    """Creates the files for the test (one file per plate.

    Args:
        sample_plates (dict): the plate for each sample.
        fam (str): the fam file.
        prefix (str): the output prefix.

    """
    logger.info("Generating plate bias files")
    samples = list(plink_utils.parse_fam(fam))

    for plate in set(sample_plates.values()):
        with open(f"{prefix}.{plate}.pheno", "w") as f:
            for sample in samples:
                sample_plate = sample_plates[(sample.fid, sample.iid)]
                status = 2 if sample_plate == plate else 1

                print(sample.fid, sample.iid, status, file=f)


def get_plates(filename: str) -> Dict[Tuple[str, str], str]:
    """Gets the plate for each samples.

    Args:
        filename (str): the file containing the plates.

    Returns:
        dict: the (FID, IID) assigned to each plate.

    """
    logger.info("Reading plate information from '%s'", filename)
    plates = {}
    with open(filename) as f:
        for line in f:
            fid, iid, plate = line.rstrip().split()
            plates[(fid, iid)] = plate
    return plates


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    if not path.isfile(args.plates):
        raise ProgramError(f"{args.plates}: no such file")


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
    group.add_argument(
        "--plates", type=str, metavar="FILE", required=True,
        help="The file containing the plate organization of each samples. "
             "Must contains three column (with no header): famID, indID and "
             "plateName.",
    )
    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--p-filter", type=float, metavar="FLOAT", default=1e-7,
        help="The significance threshold used for the plate effect. "
             "[%(default).1e]",
    )
    group.add_argument(
        "--nb-tasks", type=int, metavar="N", default=1,
        help="The number of tasks for this analysis. [%(default)d]",
    )
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="plate_bias",
        help="The prefix of the output files. [%(default)s]",
    )