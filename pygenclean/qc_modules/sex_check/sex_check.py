"""Check sample's sex using Plink."""


import argparse
import logging
import os
from os import path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from pyplink import PyPlink

from ...error import ProgramError
from ...report.summaries import SexCheckSummary
from ...utils import plink as plink_utils
from ...utils import split_extra_args, timer
from ...utils.task import execute_external_command
from ...version import pygenclean_version as __version__
from . import baf_lrr_plot, intensity_plot


SCRIPT_NAME = "sex-check"
DESCRIPTION = "Check sample's sex using Plink."
DEFAULT_OUT = "sexcheck"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> Dict[str, str]:
    """Plots the BAF and LRR of samples with sex mismatch.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the arguments as list.

    These are the steps:

    1.  Checks if there are enough markers on the chromosome ``23``.
    2.  Runs the sex check analysis using Plink.
    3.  If there are no sex problems, then quits.
    4.  Creates the recoded file for the chromosome ``23``.
    5.  Computes the heterozygosity percentage on the chromosome ``23``.
    6.  If there are enough markers on chromosome ``24`` (at least 1), creates
        the recoded file for this chromosome.
    7.  Computes the number of no call on the chromosome ``24``.
    8.  If required, plots the sex plot.
    9. If required, plots the BAF and LRR plot.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    logger.info("%s", DESCRIPTION)

    # The summary object
    summary = SexCheckSummary(args)

    # Getting the list of markers on chromosome 23
    nb_markers = len(
        plink_utils.get_markers_on_chrom(args.bfile + ".bim", {23})
    )
    if nb_markers < args.nb_chr_23:
        logger.warning(
            "There are not enough markers on chromosome 23: STOPPING NOW!",
        )
        return {
            "summary": summary,
            "usable_files": {
                "bfile": args.bfile,
            },
        }

    # Run Plink sex-check
    logger.info("Executing Plink's sex check algorithm")
    execute_external_command(
        command=[
            "plink" if args.plink_107 else "plink1.9",
            "--noweb",
            "--bfile", args.bfile,
            "--check-sex",
            "--out", args.out,
        ],
    )

    mismatches = get_sex_mismatch(
        args.out + ".sexcheck", args.female_f, args.male_f, args.out,
    )

    if len(mismatches) == 0:
        logger.info("There are no sex mismatches to investigate")

    else:
        logger.info("Computing statistics on chr23 and chr24")
        compute_statistics(
            args.bfile, mismatches, args.out,
        )

    # Intensity plot (if required)
    if args.intensity_plot:
        logger.info("Generating the intensity plot")
        intensity_plot.main(argv=[
            "--bfile", args.bfile,
            "--sex-intensities", args.sex_intensities,
            "--sex-mismatches", args.out + ".list_problem_sex",
            "--out", args.out,
        ] + args.intensity_plot_extra_args)

    if args.baf_lrr:
        logger.info("Generating LRR and BAF plot")
        dirname = args.out + ".BAF_LRR"
        if not path.isdir(dirname):
            os.mkdir(dirname)

        baf_lrr_plot.main(argv=[
            "--per-sample-baf-lrr-dir", args.per_sample_baf_lrr_dir,
            "--problematic-samples", args.out + ".list_problem_sex_ids",
            "--out", path.join(dirname, "sample"),
        ] + args.baf_lrr_extra_args)

    return {
        "summary": summary,
        "usable_files": {
            "bfile": args.bfile,
        },
    }


def compute_statistics(bfile: str, samples: Set[Tuple[str, str]],
                       prefix: str) -> None:
    """Computes the heterozygosity percentage on chromosome 23.

    Args:
        bfile (str): the prefix of the binary plink file.
        samples (set): the set of samples to check.
        prefix (str): the prefix of the output files.

    """
    with PyPlink(bfile) as bed:
        # Getting the required samples
        fam = bed.get_fam().set_index(["fid", "iid"])
        required_fam = fam.index.isin(samples)
        fam = fam.loc[required_fam, :]

        # Getting the required markers for chromosome 23
        bim = bed.get_bim()
        chr23_genotypes = get_genotypes(
            required_fam, bim.loc[bim.chrom == 23].index.tolist(), bed,
        )
        chr24_genotypes = get_genotypes(
            required_fam, bim.loc[bim.chrom == 24].index.tolist(), bed
        )

    write_heterozygosity(
        prefix + ".chr23.hetero.tsv", samples, fam, chr23_genotypes,
    )
    write_no_call(prefix + ".chr24.no_call.tsv", samples, fam, chr24_genotypes)


def write_heterozygosity(filename: str, samples: Set[Tuple[str, str]],
                         fam: pd.DataFrame, genotypes: np.ndarray) -> None:
    """Writes heterozygosity rates in output file."""
    with open(filename, "w") as f:
        print("FID", "IID", "PEDSEX", "HETERO", sep="\t", file=f)

        for i in range(len(samples)):
            # Getting the non-missing genotypes
            genos = genotypes[:, i]
            genos = genos[genos != -1]

            hetero = -9
            if genos.shape[0] > 0:
                hetero = np.sum(genos == 1) / genos.shape[0]

            # Getting the sample information
            sample = fam.iloc[i, :]

            print(*sample.name, sample.gender, hetero, sep="\t", file=f)


def write_no_call(filename: str, samples: Set[Tuple[str, str]],
                  fam: pd.DataFrame, genotypes: np.ndarray) -> None:
    """Writes the no call on chr24."""
    with open(filename, "w") as f:
        print("FID", "IID", "PEDSEX", "N_GENO", "N_MISS", "F_MISS", sep="\t",
              file=f)

        for i in range(len(samples)):
            genos = genotypes[:, i]
            nb_no_call = np.sum(genos == -1)
            pct_no_call = -9
            if genos.shape[0] > 1:
                pct_no_call = nb_no_call / genos.shape[0]

            sample = fam.iloc[i, :]

            print(*sample.name, sample.gender, genos.shape[0], nb_no_call,
                  pct_no_call, sep="\t", file=f)


def get_genotypes(samples: pd.Series, markers: List[str],
                  bed: PyPlink) -> np.ndarray:
    """Retrieves the genotype matrix."""
    # The complete genotypes
    genotypes = np.zeros(
        (len(markers), samples.sum()), dtype=np.int8,
    )

    # Parsing the data
    for i, marker in enumerate(markers):
        marker_genotypes = bed.get_geno_marker(marker)[samples]
        genotypes[i, :] = marker_genotypes

    return genotypes


def get_sex_mismatch(filename: str, female_f: float, male_f: float,
                     prefix: str) -> Set[Tuple[str, str]]:
    """Gets the sex mismatches.

    Args:
        filename (str): the name of Plink's sex check file.
        female_f (float): the F threshold for females.
        male_f (float): the F threshold for males.
        prefix (str): the output prefix.

    Returns:
        set: the set of samples with sex mismatch (if any).

    """
    # pylint: disable=no-member

    if not path.isfile(filename):
        raise ProgramError(f"{filename}: something went wrong with PLink")

    df = pd.read_csv(filename, sep=r"\s+")

    # Getting some statistics
    #   - the problems (i.e. != OK)
    #   - males (according to PEDSEX, i.e. == 1)
    #   - females (according to PEDSEX, i.e. == 2)
    #   - unknown PEDSEX (i.e. == 0)
    #   - unknown SNPSEX (i.e. == 0)
    probs = df.STATUS != "OK"
    males = df.PEDSEX == 1
    females = df.PEDSEX == 2
    pedsex_unknown = df.PEDSEX == 0
    snpsex_unknown = df.SNPSEX == 0

    # Finding the issues
    # The mismatches are the samples where PEDSEX != SNPSEX
    mismatches = df.PEDSEX != df.SNPSEX

    # We want to rescue females where SNPSEX is unknown and lower than F
    ok_females = mismatches & females & snpsex_unknown & (df.F < female_f)

    # We want to rescue males where SNPSEX is unknown and higher than F
    ok_males = mismatches & males & snpsex_unknown & (df.F > male_f)

    logger.info("sex-check summary")
    logger.info("  - %d total problems", probs.sum())
    logger.info("  - %d pedsex unknown", pedsex_unknown.sum())
    logger.info("  - %d female F < %.1f", ok_females.sum(), female_f)
    logger.info("  - %d male F > %.1f", ok_males.sum(), male_f)

    # Keeping only the problematic samples
    df = df.loc[mismatches & (~ (ok_females | ok_males)), :]
    logger.info("  - %d mismatches kept", df.shape[0])

    # Saving the list of problems
    df.to_csv(prefix + ".list_problem_sex", index=False, sep="\t")

    # Saving the list of IDs with problems
    df.loc[:, ["FID", "IID"]].to_csv(
        prefix + ".list_problem_sex_ids", index=False, header=False, sep="\t",
    )

    return set(
        tuple(values) for _, values in df.loc[:, ["FID", "IID"]].iterrows()
    )


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options to check.

    """
    # Check if we have the tped and the tfam files
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    # Ceck the number of markers on chromosome 23
    if args.nb_chr_23 < 0:
        raise ProgramError(
            f"{args.nb_chr_23}: number of markers on chr23 must be "
            f"positive",
        )

    # Intensity plot
    if args.intensity_plot:
        if args.sex_intensities is None:
            raise ProgramError(
                "Asking for intensity plot, but no '--intensities' provided"
            )

        args.intensity_plot_extra_args = split_extra_args(
            args.intensity_plot_extra_args,
        )

    # BAF and LRR plot
    if args.baf_lrr:
        if args.per_sample_baf_lrr_dir is None:
            raise ProgramError(
                "Asking for BAF & LRR plot, but no '--intensity-dir' provided"
            )

        args.baf_lrr_extra_args = split_extra_args(args.baf_lrr_extra_args)


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
        help="The input file prefix (will find the Plink binary files by "
             "appending the prefix to the .bed, .bim, and .fam files, "
             "respectively.",
    )
    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--female-f", type=float, metavar="FLOAT", default=0.3,
        help="The female F threshold. [< %(default)f]",
    )
    group.add_argument(
        "--male-f", type=float, metavar="FLOAT", default=0.7,
        help="The male F threshold. [> %(default)f]",
    )
    group.add_argument(
        "--nb-chr-23", type=int, metavar="INT", default=50,
        help="The minimum number of markers on chromosome 23 before computing "
             "Plink's sex check [%(default)d]",
    )
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )

    group = parser.add_argument_group("Intensity plot")
    group.add_argument(
        "--intensity-plot", action="store_true",
        help="Create the sex plot (summarized chr Y intensities in "
             "function of summarized chr X intensities) for problematic "
             "samples.",
    )
    group.add_argument(
        "--sex-intensities", type=str, metavar="FILE",
        help="A file containing allele intensities for each of the markers "
             "located on the X and Y chromosome for the sex plot. Note "
             "that all samples need to be in this file.",
    )
    group.add_argument(
        "--intensity-plot-format", type=str, metavar="FORMAT",
        default="png", choices=["png", "ps", "pdf"],
        help="The output file format for the sex plot (png, ps or pdf "
             "formats are available). [default: %(default)s]",
    )
    group.add_argument(
        "--intensity-plot-extra-args", type=str, metavar="EXTRA",
        help="Extra arguments for intensity-plot (see help for more "
             "information).",
    )

    group = parser.add_argument_group("LRR and BAF plot")
    group.add_argument(
        "--baf-lrr", action="store_true",
        help="Create the LRR and BAF plot for problematic samples",
    )
    group.add_argument(
        "--per-sample-baf-lrr-dir", type=str, metavar="DIR",
        help="Directory containing BAF and LRR values for every samples (one "
             "file per sample).",
    )
    group.add_argument(
        "--baf-lrr-format", type=str, metavar="FORMAT",
        default="png", choices=["png", "ps", "pdf"],
        help="The output file format for the LRR and BAF plot (png, ps or pdf "
             " formats are available). [%(default)s]",
    )
    group.add_argument(
        "--baf-lrr-extra-args", type=str, metavar="EXTRA",
        help="Extra arguments for baf-lrr-plot (see help for more "
             "information).",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output files")
    group.add_argument(
        "--out", type=str, metavar="FILE", default=DEFAULT_OUT,
        help="The prefix of the output files (which will be a Plink binary "
             "file. [%(default)s]",
    )
