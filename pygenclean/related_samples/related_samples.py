"""Finds related samples according to IBS (if any)."""


import logging
import argparse

from ..error import ProgramError

from ..utils import task
from ..utils import plink as plink_utils

from ..version import pygenclean_version as __version__


SCRIPT_NAME = "related-samples"
DESCRIPTION = "Finds related samples according to IBS values."

_IBS2_RATIO_DEFAULT = 0.8
_INDEP_PAIRWISE_R2_DEFAULT = "0.1"

logger = logging.getLogger(__name__)


def main(args=None, argv=None):
    """Finds related samples according to IBS (if any).

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1. Uses Plink to extract markers according to LD.
    2. Checks if there is enough markers after pruning. If not, then quits.
    3. Extract markers according to LD.
    4. Runs Plink with the ``genome`` option. Quits here if the user asker only
       for the ``genome`` file.
    5. Finds related individuals and gets values for plotting.
    6. Plots ``Z1`` in function of ``IBS2 ratio`` for related individuals.
    7. Plots ``Z2`` in function of ``IBS2 ratio`` for related individuals.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    markers_to_extract = select_markers_by_ld(
        args.bfile, args.maf, args.indep_pairwise, args.out,
    )

    # Counting the number of markers
    with open(markers_to_extract) as f:
        nb_markers = len(f.read().splitlines())
    if nb_markers < args.min_nb_snp:
        logger.warning("Only %d markers on autosome, stopping", nb_markers)
        return

    # Extracting the markers
    plink_utils.extract_markers(
        bfile=args.bfile, extract=markers_to_extract,
        out=args.out + ".pruned_data",
    )


def select_markers_by_ld(bfile, maf_threshold, indep_pairwise, out):
    """Selects markers according to LD.

    Args:
        bfile (str): the prefix of the plink files.
        maf_threshold (str): the MAF threshold.
        indep_pairwise (list): the values for indep pairwise selection.
        out (str): the output file prefix.

    Returns:
        str: the name of the file containing the markers to extract.

    Note:
        This script uses Plink version 2.

    """
    logger.info("Extract markers according to LD")

    out_prefix = out + ".pruning_" + indep_pairwise[2]

    command = [
        "plink2",
        "--bfile", bfile,
        "--maf", maf_threshold,
        "--out", out_prefix,
        "--indep-pairwise",
    ] + indep_pairwise
    task.execute_external_command(command)

    # Finding autosomal markers
    autosomes = set(range(1, 23))
    autosomal_markers = plink_utils.get_markers_on_chrom(
        bfile + ".bim", autosomes,
    )

    # Reading the pruned markers
    with open(out_prefix + ".prune.in", "r") as f:
        pruned_markers = set(f.read().splitlines())

    # Writing the pruned markers on autosomes
    with open(out_prefix + ".prune.in.autosomal", "w") as f:
        print(*autosomal_markers & pruned_markers, sep="\n", file=f)

    return out_prefix + ".prune.in.autosomal"


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: missing plink files")


def parse_args(argv=None):
    """Parses the arguments and function."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments
    add_args(parser)

    if argv is None:
        return parser.parse_args()
    return parser.parse_args(argv)


def add_args(parser):
    """Adds argument to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input files")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively.)",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--genome-only", action="store_true",
        help="Only create the genome file",
    )
    group.add_argument(
        "--min-nb-snp", type=int, metavar="INT", default=10000,
        help="The minimum number of markers needed to compute IBS values. "
             "[%(default)d]",
    )
    group.add_argument(
        "--indep-pairwise", type=str, metavar="STR", nargs=3,
        default=["50", "5", _INDEP_PAIRWISE_R2_DEFAULT],
        help="Three numbers: window size, window shift and the r2 threshold. "
             "[%(default)s]",
    )
    group.add_argument(
        "--maf", type=str, metavar="FLOAT", default="0.05",
        help="Restrict to SNPs with MAF >= threshold. [%(default)s]",
    )
    group.add_argument(
        "--ibs2-ratio", type=float, metavar="FLOAT",
        default=_IBS2_RATIO_DEFAULT,
        help="The initial IBS2* ratio (the minimum value to show in the "
             "plot. [%(default).1f]",
    )
    group.add_argument(
        "--nb-tasks", type=int, metavar="N", default=1,
        help="The number of tasks for this analysis. [%(default)d]",
    )
    group.add_argument(
        "--lines-per-task", type=int, metavar="INT", default=100,
        help="The number of lines per task. [%(default)d]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="ibs",
        help="The prefix of the output files. [%(default)s]",
    )
