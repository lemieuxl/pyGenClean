"""Finds related samples according to IBS (if any)."""


import argparse
import logging
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ...error import ProgramError
from ...report.summaries import RelatedSamplesSummary
from ...utils import count_lines
from ...utils import plink as plink_utils
from ...utils import task, timer
from ...version import pygenclean_version as __version__
from . import merge_related_samples


SCRIPT_NAME = "related-samples"
DESCRIPTION = "Finds related samples according to IBS values."
DEFAULT_OUT = "ibs"


_IBS2_RATIO_DEFAULT = 0.8
_INDEP_PAIRWISE_R2_DEFAULT = "0.1"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> Dict[str, str]:
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

    summary = RelatedSamplesSummary(args)

    logger.info("%s", DESCRIPTION)

    markers_to_extract = select_markers_by_ld(
        args.bfile, args.maf, args.indep_pairwise, args.out, args.plink_107,
    )

    # Counting the number of markers
    nb_markers = count_lines(markers_to_extract)

    if nb_markers < args.min_nb_snp:
        logger.warning("Only %d markers on autosome, stopping", nb_markers)
        return {
            "summary": summary,
            "usable_files": {
                "bfile": args.bfile,
            },
        }

    # Extracting the markers
    plink_utils.subset(
        bfile=args.bfile,
        out=get_prefix_for_genome(args.out),
        markers=markers_to_extract,
        marker_subset_type="extract",
        use_original_plink=args.plink_107,
    )

    # Creating the genome file
    create_genome(bfile=get_prefix_for_genome(args.out), args=args)

    # We only need the genome file to be created
    if args.genome_only:
        return {
            "summary": summary,
            "usable_files": {
                "bfile": args.bfile,
            },
        }

    # Getting the related samples
    related = get_related_samples(args.out + ".genome.genome.gz", args.out,
                                  args.ibs2_ratio)

    if related is None:
        logger.info("There are no related samples in the dataset")
        return {
            "summary": summary,
            "usable_files": {
                "bfile": args.bfile,
            },
        }

    # Merge related samples
    merge_related_samples.main(
        argv=[
            "--ibs-related", args.out + ".related_individuals",
            "--out", args.out,
        ],
    )

    # Plotting the related samples (z1 and z2)
    plot_related_samples(related, args.out, args)

    return {
        "summary": summary,
        "usable_files": {
            "bfile": args.bfile,
            "discarded": args.out + ".discarded_related_individuals",
        },
    }


def get_prefix_for_genome(prefix: str) -> str:
    """Generate the prefix for the creation of the genome file"""
    return prefix + ".pruned_data"


def create_genome(bfile: str, args: argparse.Namespace) -> None:
    """Runs the genome command from plink.

    Args:
        bfile (str): the input file prefix.
        options (argparse.Namespace): the options.

    Returns:
        str: the name of the ``genome`` file.

    Runs Plink with the ``genome`` option. If the user asks for SGE
    (``options.sge`` is True), a frequency file is first created by plink.
    Then, the input files are split in ``options.line_per_file_for_sge`` and
    Plink is called (using the ``genome`` option) on the cluster using SGE
    (:py:func:`runGenomeSGE`). After the analysis, Plink's output files and
    logs are merged using :py:func:`mergeGenomeLogFiles`.

    """
    logger.info("Creating the genome file")

    # Generating the genome files
    command = [
        "plink" if args.plink_107 else "plink1.9",
        "--noweb",
        "--bfile", bfile,
        "--out", args.out + ".genome",
    ]

    if args.plink_107:
        command.extend([
            "--Z-genome",
            "--genome-full",
        ])

    else:
        command.extend([
            "--genome", "gz", "full",
            "--threads", str(args.nb_threads),
        ])

    task.execute_external_command(command)


def get_related_samples(
    filename: str,
    out: str,
    ibs2_ratio_threshold: float,
) -> Optional[pd.DataFrame]:
    """Get the related samples according to IBS2* ratio."""
    logger.info("Reading genome file to get related samples")

    # This file is pretty big, so we neet to parse it chuinks by chunks
    reader = pd.read_csv(
        filename,
        delim_whitespace=True,
        dtype={"FID1": str, "IID1": str, "FID2": str, "IID2": str},
        chunksize=10e6,
    )

    # This list we be concatenated into a single DataFrame
    final: pd.DataFrame = []

    for df in reader:
        df["IBS2_ratio"] = df.HETHET / (df.HOMHOM + df.HETHET)

        # If IBS2* ratio is bigger than the threshold, they might be related
        df = df.loc[df.IBS2_ratio > ibs2_ratio_threshold, :]

        # Are there any related samples?
        if df.shape[0] == 0:
            continue

        df["status"] = "unknown"
        df["code"] = 5

        # Full sibs
        subset = (
            ((df.Z0 >= 0.17) & (df.Z0 <= 0.33))
            & ((df.Z1 >= 0.40) & (df.Z1 <= 0.60))
        )
        df.loc[subset, "status"] = "full-sibs"
        df.loc[subset, "code"] = 1

        # Half sibs, grand-parent child, uncle newphew
        subset = (
            (df.Z0 >= 0.4) & (df.Z0 <= 0.6) & (df.Z1 >= 0.4) & (df.Z1 <= 0.6)
        )
        df.loc[subset, "status"] = "half-sibs;grand-parent-child;uncle-nephew"
        df.loc[subset, "code"] = 2

        # Parent child
        subset = (df.Z0 <= 0.05) & (df.Z1 >= 0.95) & (df.Z2 <= 0.05)
        df.loc[subset, "status"] = "parent-child"
        df.loc[subset, "code"] = 3

        # Twins
        subset = (df.Z0 <= 0.05) & (df.Z1 <= 0.05) & (df.Z2 >= 0.95)
        df.loc[subset, "status"] = "twins"
        df.loc[subset, "code"] = 4

        final.append(df)

    # Do we have some related samples?
    if len(final) == 0:
        return None

    # Creating the final data frame
    final = pd.concat(final, ignore_index=True)

    # Saving to file
    final.to_csv(out + ".related_individuals", index=False, sep="\t",
                 na_rep="NA")

    return final


def select_markers_by_ld(
    bfile: str,
    maf_threshold: str,
    indep_pairwise: List[str],
    out: str,
    use_original_plink: bool = False,
) -> str:
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
        "plink" if use_original_plink else "plink1.9",
        "--noweb",
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


def plot_related_samples(
    data: pd.DataFrame,
    out: str,
    args: argparse.Namespace,
) -> None:
    """Plot Z1 and Z2 in function of IBS2* ratio."""
    # Creating the plots
    for z_col, y_label in zip(("Z1", "Z2"), ("$Z_1$", "$Z_2$")):
        figure, axe = plt.subplots(1, 1)

        # Annotation
        axe.set_title(
            f"{data.shape[0]} pairs with "
            r"$IBS2^\ast_{ratio} >"
            f"{args.ibs2_ratio}$"
        )
        axe.set_xlabel(r"$IBS2^\ast_{ratio}$")
        axe.set_ylabel(y_label)

        # Plotting the data (there are 5 error codes)
        # Code 1 (full sibs)
        code_1 = data.code == 1
        axe.scatter(
            x=data.loc[code_1, "IBS2_ratio"],
            y=data.loc[code_1, z_col],
            color="#CC0000",
            s=5,
            label=f"Full sibs (n={np.count_nonzero(code_1)})",
            zorder=2,
        )

        # Code 2 (half sibs)
        code_2 = data.code == 2
        axe.scatter(
            x=data.loc[code_2, "IBS2_ratio"],
            y=data.loc[code_2, z_col],
            color="#0099CC",
            s=5,
            label=f"Half sibs, grand-parent-child or uncle-nephew "
                  f"(n={np.count_nonzero(code_2)})",
            zorder=2,
        )

        # Code 3 (parent-child)
        code_3 = data.code == 3
        axe.scatter(
            x=data.loc[code_3, "IBS2_ratio"],
            y=data.loc[code_3, z_col],
            color="#FF8800",
            s=5,
            label=f"Parent-child (n={np.count_nonzero(code_3)})",
            zorder=2,
        )

        # Code 4 (twins / duplicated samples)
        code_4 = data.code == 4
        axe.scatter(
            x=data.loc[code_4, "IBS2_ratio"],
            y=data.loc[code_4, z_col],
            color="#9933CC",
            s=5,
            label=f"Twins or duplicated samples "
                  f"(n={np.count_nonzero(code_4)})",
            zorder=2,
        )

        # Code 5 (unknown)
        code_5 = data.code == 5
        axe.scatter(
            x=data.loc[code_5, "IBS2_ratio"],
            y=data.loc[code_5, z_col],
            color="#669900",
            s=5,
            label=f"Unknown (n={np.count_nonzero(code_5)})",
            zorder=1,
        )

        # Setting the X and Y limits
        axe.set_xlim((args.ibs2_ratio - 0.01, 1.01))
        axe.set_ylim((-0.01, 1.01))

        # Adding the legend
        axe.legend(fontsize=8)

        # Modifying the spines
        axe.xaxis.set_ticks_position("bottom")
        axe.yaxis.set_ticks_position("left")
        axe.spines["top"].set_visible(False)
        axe.spines["right"].set_visible(False)

        # Saving the figure
        plt.savefig(out + f".related_individuals_{z_col.lower()}.png")
        plt.close(figure)


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")


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
        "--nb-threads", type=int, metavar="N", default=1,
        help="The number of threads for this analysis (no effect when using "
             "plink 1.07). [%(default)d]",
    )
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07). Note that this will be slow, "
             "as Plink 1.07 doesn't support multi threading.",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default=DEFAULT_OUT,
        help="The prefix of the output files. [%(default)s]",
    )
