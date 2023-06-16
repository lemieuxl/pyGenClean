"""Check for sample contamination using BAFRegress."""


import argparse
import logging
import math
from os import path
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

from ...error import ProgramError
from ...report.summaries import ContaminationSummary
from ...utils import plink as plink_utils
from ...utils import task as task_utils
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "contamination"
DESCRIPTION = "Check for sample contamination using BAFRegress."
DEFAULT_OUT = "contamination"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> Dict[str, str]:
    """Check for sample contamination using BAFRegress.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1. Prints the options.
    2. Compute frequency using Plink.
    3. Runs bafRegress.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    logger.info("%s", DESCRIPTION)

    # Checks the sample raw data
    logger.info("Checking the raw data files")
    sample_files = check_sample_files(args.bfile + ".fam", args.raw_dir)

    # Finds the markers to extract
    logger.info("Creating extraction list (autosome only)")
    create_extraction_file(args.bfile + ".bim", args.out)

    # Run plink
    logger.info("Computing frequency using Plink")
    plink_utils.compute_freq(
        bfile=args.bfile,
        out=args.out,
        extract=args.out + ".to_extract",
        use_original_plink=args.plink_107,
    )

    run_bafregress(
        filenames=sample_files,
        out=args.out,
        extract=args.out + ".to_extract",
        freq=args.out + ".frq",
        args=args,
    )

    # Writing contaminated samples
    write_contaminated_samples(
        filename=args.out + ".bafRegress",
        out=args.out,
        threshold=args.estimate_threshold,
    )

    return {
        "summary": ContaminationSummary(args),
        "usable_files": {
            "bfile": args.bfile,
            "contaminated": args.out + ".contaminated_samples",
        },
    }


def check_sample_files(fam: str, raw_dirname: str) -> Set[str]:
    """Checks the raw sample files.

    Args:
        fam (str): the name of the FAM file.
        raw_dirname (str): the name of the directory containing the raw file.

    Returns:
        set: the set of all the sample files that are compatible with the FAM
             file.

    """
    # Reading the sample identification number from the FAM file
    fam_samples = set(row.iid for row in plink_utils.parse_fam(fam))

    # Checking the files in the raw directory
    sample_files = set()
    all_samples = set()
    for filename in Path(raw_dirname).glob("*"):
        # Keeping only the files
        if not filename.is_file():
            continue

        sample = path.splitext(filename.name)[0]
        all_samples.add(sample)

        if sample not in fam_samples:
            logger.warning("%s: sample not in FAM file", sample)
        else:
            sample_files.add(str(filename))

    for sample in fam_samples - all_samples:
        logger.warning("%s: sample not in raw directory", sample)

    if len(sample_files) == 0:
        raise ProgramError("no sample left for analysis")

    return sample_files


def create_extraction_file(bim: str, out: str):
    """Creates an extraction file (keeping only markers on autosomes).

    Args:
        bim (str): the name of the BIM file.
        out (str): the prefix for the output file.

    """
    # Fetching the markers on the autosomes
    autosomal = plink_utils.get_markers_on_chrom(bim, set(range(1, 23)))

    with open(out + ".to_extract", "w") as f:
        print(*autosomal, sep="\n", file=f)

    logger.info("  - %s markers will be used for contamination estimation",
                f"{len(autosomal):,d}")


def run_bafregress(filenames: Set[str], out: str, extract: str, freq: str,
                   args: argparse.Namespace) -> None:
    """Runs the bafRegress function.

    Args:
        filenames (set): the set of all sample files.
        out (str): the output prefix.
        extract (str): the name of the markers to extract.
        freq (str): the name of the file containing the frequency.
        args (argparse.Namespace): the other options.

    """
    # The base command
    base_command = [
        "bafRegress.py",
        "estimate",
        "--freqfile", freq,
        "--freqcol", "2,5",
        "--extract", extract,
        "--colsample", args.colsample,
        "--colmarker", args.colmarker,
        "--colbaf", args.colbaf,
        "--colab1", args.colab1,
        "--colab2", args.colab2,
    ]

    filenames = list(filenames)

    # Creating equally numbered chunks (whith maximum number of samples)
    results = []
    for j in range(0, len(filenames), args.max_samples_per_job * args.nb_cpu):
        sample_chunk = filenames[j: j + args.max_samples_per_job * args.nb_cpu]

        nb_samples_per_chunk = math.ceil(len(sample_chunk) / args.nb_cpu)
        chunks = [
            sample_chunk[i:i + nb_samples_per_chunk]
            for i in range(0, len(sample_chunk), nb_samples_per_chunk)
        ]

        # Creating and executing the commands
        results += task_utils.execute_external_commands(
            [base_command + chunk for chunk in chunks], args.nb_cpu,
        )

    # Saving the output
    header = None
    with open(out + ".bafRegress", "w") as f:
        for result in results:
            content = result.splitlines()
            if header is None:
                header = content[0]
                print(*content, sep="\n", file=f)

            else:
                if content[0] != header:
                    raise ProgramError("different header from BAFregress")
                print(*content[1:], sep="\n", file=f)


def write_contaminated_samples(filename: str, out: str,
                               threshold: float) -> None:
    """Write contaminated samples to file."""
    # Reading the contaminated samples
    df = pd.read_csv(filename, sep="\t")
    contaminated = df.loc[df.estimate > threshold, "sample"]

    # Printing to file
    with open(out + ".contaminated_samples", "w") as f:
        for sample in contaminated.to_list():
            print(sample, sample, file=f)


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    # Checking that the raw directory exists
    if not path.isdir(args.raw_dir):
        raise ProgramError(f"{args.raw_dir}: no such directory")

    # Checking the estimate is between 0 and 1
    if not 0 <= args.estimate_threshold <= 1:
        raise ProgramError(f"{args.estimate_threshold}: not between 0 and 1")


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}"
    )

    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Adds argument to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively).",
    )

    # The raw data options
    group = parser.add_argument_group("Raw Data")
    group.add_argument(
        "--raw-dir", type=str, metavar="DIR", required=True,
        help="Directory containing the raw data (one file per sample, where "
             "the name of the file (minus the extension) is the sample "
             "identification number.",
    )
    group.add_argument(
        "--colsample", type=str, metavar="COL", default="Sample Name",
        help="The sample column. [default: %(default)s]",
    )
    group.add_argument(
        "--colmarker", type=str, metavar="COL", default="SNP Name",
        help="The marker column. [default: %(default)s]",
    )
    group.add_argument(
        "--colbaf", type=str, metavar="COL", default="B Allele Freq",
        help="The B allele frequency column. [default: %(default)s]",
    )
    group.add_argument(
        "--colab1", type=str, metavar="COL", default="Allele1 - AB",
        help="The AB Allele 1 column. [default: %(default)s]",
    )
    group.add_argument(
        "--colab2", type=str, metavar="COL", default="Allele2 - AB",
        help="The AB Allele 2 column. [default: %(default)s]",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--estimate-threshold", type=float, metavar="FLOAT", default=0.01,
        help="The estimate threshold for which a sample is considered "
             "contaminated. [>%(default).2f]",
    )
    group.add_argument(
        "--nb-cpu", type=int, metavar="NB", default=1,
        help="The number of CPU to use. [default: %(default)s]",
    )
    group.add_argument(
        "--max-samples-per-job", type=int, metavar="NB", default=100,
        help="The maximum number of samples per task. [%(default)d]",
    )
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07). Note that this will be slow, "
             "as Plink 1.07 doesn't support multi threading.",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default=DEFAULT_OUT,
        help="The prefix of the output files. [default: %(default)s]",
    )
