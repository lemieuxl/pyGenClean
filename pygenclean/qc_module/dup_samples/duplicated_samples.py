"""Verify completion and concordance of duplicated samples."""


import argparse
import collections
import logging
import random
from os import path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from pyplink import PyPlink

from ...error import ProgramError
from ...utils import plink as plink_utils
from ...utils import timer
from ...version import pygenclean_version as __version__


SCRIPT_NAME = "duplicated-samples"
DESCRIPTION = "Verify completion and concordance of duplicated samples."


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Check for duplicated samples.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as a list.

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    # Find duplicated samples
    logger.info("Reading the duplicated samples")
    duplicated_samples = read_duplicates(args.duplicated_samples)

    # Generaeting the files contianing only the unique samples
    with open(args.out + ".duplicated_samples", "w") as f:
        for duplicates in duplicated_samples.values():
            for sample in duplicates:
                print(*sample, file=f)

    # Generating the files with the unique samples
    logger.info("Generating the plink files with unique samples")
    plink_utils.subset_samples(
        bfile=args.bfile,
        samples=args.out + ".duplicated_samples",
        out=args.out + ".unique_samples",
        subset_type="remove",
        use_original_plink=args.plink_107,
    )

    # Generating the files with the duplicated samples
    logger.info("Generating the plink files with duplicated samples")
    plink_utils.subset_samples(
        bfile=args.bfile,
        samples=args.out + ".duplicated_samples",
        out=args.out + ".duplicated_samples",
        subset_type="keep",
        use_original_plink=args.plink_107,
    )

    # Compute statistics
    logger.info(
        "Computing the completion and concordance of duplicated samples",
    )
    completion, concordance, sample_order = compute_statistics(
        bfile=args.out + ".duplicated_samples",
        dup_samples=duplicated_samples,
        out=args.out,
    )

    # Choose the best duplicates
    logger.info("Choosing the best duplicates")
    choose_best_duplicates(
        dup_samples=sample_order, completion=completion,
        concordance=concordance, out=args.out,
    )

    # Generating the files with the chosen samples
    logger.info("Generating the plink files with only the chosen samples")
    plink_utils.subset_samples(
        bfile=args.out + ".duplicated_samples",
        samples=args.out + ".chosen_samples.info",
        out=args.out + ".chosen_samples",
        subset_type="keep",
        use_original_plink=args.plink_107,
    )

    # Merging the two files (chosen and unique samples)
    plink_utils.merge_files(
        prefixes=[args.out + ".unique_samples", args.out + ".chosen_samples"],
        out=args.out + ".final",
        use_original_plink=args.plink_107
    )


def choose_best_duplicates(
    dup_samples: Dict[str, List[Tuple[str, str]]],
    *,
    completion: Dict[str, np.ndarray],
    concordance: Dict[str, np.ndarray],
    out: str,
) -> None:
    """Choose the best duplicates according to the completion rate.

    Args:
        dup_samples (dict): the updated position of the samples in the tped
                            containing only duplicated samples.
        completion (dict): the completion of each of the duplicated samples.
        concordance (dict): the concordance of every duplicated samples.
        out (str): the prefix of all the files.

    Returns:
        tuple: a tuple where the first element is a dictionary of the chosen
               samples and the second one is the mean concordance.

    These are the steps to find the best duplicated sample:

    1. Sort the list of concordances.
    2. Sort the list of completions.
    3. Choose the best of the concordance and put in a set.
    4. Choose the best of the completion and put it in a set.
    5. Compute the intersection of the two sets. If there is one sample or
       more, then randomly choose one sample.
    6. If the intersection doesn't contain at least one sample, redo steps 3
       and 4, but increase the number of chosen best by one. Redo step 5 and 6
       (if required).

    The chosen samples are written in ``prefix.chosen_samples.info``. The rest
    are written in ``prefix.excluded_samples.info``.

    """
    # For each duplicated sample
    chosen = {}
    for sample, duplicates in dup_samples.items():
        # Getting the completion for those duplicated samples and sorting it
        curr_completion = completion[sample]
        sorted_completion_idx = np.argsort(curr_completion)

        # Computing the mean concordance for each sample and sorting it
        curr_concordance = (
            (concordance[sample].sum(axis=1) - np.diag(concordance[sample]))
            / (concordance[sample].shape[0] - 1)
        )
        sorted_concordance_idx = np.argsort(curr_concordance)

        # Trying to find the best duplicate to keep
        nb_to_check = 1
        chosen_index = None
        while nb_to_check <= len(duplicates):
            # Getting the `nb_to_check` best value (higher to lower)
            completion_value = curr_completion[
                sorted_completion_idx[nb_to_check*-1]
            ]
            concordance_value = curr_concordance[
                sorted_concordance_idx[nb_to_check*-1]
            ]

            # Getting the indexes to consider
            completion_to_consider = set(
                np.where(curr_completion >= completion_value)[0]
            )
            concordance_to_consider = set(
                np.where(curr_concordance >= concordance_value)[0]
            )

            # Getting the intersection of the indexes
            to_consider = concordance_to_consider & completion_to_consider
            if len(to_consider) >= 1:
                chosen_index = random.choice(list(to_consider))
                break
            nb_to_check += 1

        if chosen_index is None:
            raise ProgramError(
                f"Could not choose the best sample ID for {sample}",
            )

        chosen[sample] = duplicates[chosen_index]

    # Printing the chose and excluded samples
    with open(out + ".chosen_samples.info", "w") as f:
        for fid, iid in chosen.values():
            logger.info("  - %s/%s was chosen", fid, iid)
            print(fid, iid, file=f)

    # Printing the excluded samples
    with open(out + ".excluded_samples.info", "w") as f:
        for sample, duplicates in dup_samples.items():
            for fid, iid in duplicates:
                if (fid, iid) != chosen[sample]:
                    print(fid, iid, file=f)


def print_concordance(concordance: Dict[str, np.ndarray], out: str) -> None:
    """Print the concordance.

    Args:
        concordance (dict): the concordance of each sample.
        out (str): the prefix of all the files.

    The concordance is the number of genotypes that are equal when comparing a
    duplicated samples with another one, divided by the total number of
    genotypes (excluding genotypes that are no call [*i.e.* ``0``]). If a
    duplicated sample has 100% of no calls, the concordance will be zero.

    The file ``prefix.concordance`` will contain :math:`N \\times N` matrices
    for each set of duplicated samples.

    """
    with open(out + ".concordance", "w") as f:
        for sample, sample_concordance in concordance.items():
            print("#" + sample, file=f)
            for i in range(sample_concordance.shape[0]):
                print(*sample_concordance[i], sep="\t", file=f)


def print_statistics(
    dup_samples: Dict[str, np.ndarray],
    *,
    completion: Dict[str, np.ndarray],
    concordance: Dict[str, np.ndarray],
    out: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Print the statistics in a file.

    Args:
        completion (dict): the completion of each duplicated samples.
        concordance (dict): the concordance of each duplicated samples.
        tpedSamples: the updated position of the samples in the tped
                     containing only duplicated samples.
        oldSamples: the original duplicated sample positions.
        out (str): the prefix of the output files.

    Returns
        dict: the concordance for each duplicated samples.

    Prints the statistics (completion of each samples and pairwise concordance
    between duplicated samples) in a file (``prefix.summary``).

    """
    # The percentages
    sample_concordance = {
        sample: sample_concordance[0] / sample_concordance[1]
        for sample, sample_concordance in concordance.items()
    }
    sample_completion = {
        sample: sample_completion[0] / sample_completion[1]
        for sample, sample_completion in completion.items()
    }

    with open(out + ".summary", "w") as f:
        print("sample", "fid", "iid", "% completion", "completion",
              "mean concordance", sep="\t", file=f)

        for sample, duplicates in dup_samples.items():
            for i, (fid, iid) in enumerate(duplicates):
                # The completion
                dup_completion = completion[sample][:, i]

                # The concordance
                dup_concordance = sample_concordance[sample]

                # Printing for each duplicates
                print(sample, fid, iid, sample_completion[sample][i],
                      f"{dup_completion[0]}/{dup_completion[1]}",
                      np.mean(np.delete(dup_concordance[i], i)),
                      sep="\t", file=f)

    return sample_concordance, sample_completion


def compute_statistics(
    bfile: str,
    dup_samples: Dict[str, Set[Tuple[str, str]]],
    out: str,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray],
           Dict[str, List[Tuple[str, str]]]]:
    """Computes the completion and concordance of each samples.

    Args:
        bfile (str): the prefix of the file containing duplicated samples.
        samples (dict): the duplicated samples.
        out (str): the output prefix.

    Returns:
        tuple: a tuple containing the completion (:py:class:`numpy.array`) as
               first element, and the concordance (:py:class:`dict`) as the
               second element and the sample order in the first element.

    """
    # The completion
    completion = {
        sample: np.zeros((2, len(duplicates)), dtype=int)
        for sample, duplicates in dup_samples.items()
    }

    # The concordance
    concordance = {
        sample: np.zeros((2, len(duplicates), len(duplicates)), dtype=int)
        for sample, duplicates in dup_samples.items()
    }

    # The sample mask and order of appearance
    sample_masks = {}
    sample_order = {}

    differences = []
    with PyPlink(bfile) as bed:
        bim = bed.get_bim()
        fam = bed.get_fam().set_index(["fid", "iid"], verify_integrity=True)

        # Getting the samples' information  and mask
        sexes = {}
        for sample, duplicates in dup_samples.items():
            # Saving the sample mask and order of appearance
            sample_masks[sample] = fam.index.isin(duplicates)
            sample_order[sample] = fam.index[sample_masks[sample]].to_list()

            # Getting the samples' sex
            sample_sex = fam.loc[sample_masks[sample], "gender"].unique()
            if sample_sex.shape[0] != 1:
                raise ProgramError(f"{sample}: multiple sexes: {sample_sex}")
            sexes[sample] = sample_sex[0]

        for marker, genotypes in bed.iter_geno():
            # The marker's chromosome
            chromosome = bim.loc[marker, "chrom"]

            # The marker's possible genotypes

            for sample, sample_mask in sample_masks.items():
                # We skip if we have a female on the Y chromosome
                if chromosome == 24 and sexes[sample] == 2:
                    continue

                # The samples' genotypes
                sample_genotypes = genotypes[sample_mask]

                # Not-missing genotypes
                not_missing = sample_genotypes != -1

                # Computing the completion
                mask = not_missing
                if sexes[sample] == 1 and chromosome in (23, 24):
                    # Males on sexual chromosome should be homozygous
                    mask = mask & (sample_genotypes != 1)
                completion[sample][0][mask] += 1
                completion[sample][1] += 1

                # Computing the concordance
                identity = sample_genotypes.reshape(-1, 1) == sample_genotypes
                mask = not_missing.reshape(-1, 1) & not_missing
                concordance[sample][0][identity & mask] += 1
                concordance[sample][1][mask] += 1

                # The unique samples' genotypes (removing missing)
                unique_genotypes = np.unique(sample_genotypes[not_missing])

                # Is there more than one unique genotype (i.e. difference)
                if unique_genotypes.shape[0] > 1:
                    # Getting the ACGT genotype map
                    acgt_geno_map = get_acgt_geno_map(
                        *bim.loc[marker, ["a1", "a2"]],
                    )

                    for fid_iid, geno in zip(fam.index[sample_mask],
                                             genotypes[sample_mask]):
                        differences.append((
                            marker, sample, *fid_iid, acgt_geno_map[geno],
                        ))

    # Printing the differences
    with open(out + ".diff", "w") as f:
        print("name", "sample", "fid", "iid", "genotype", sep="\t",
              file=f)
        for difference in differences:
            print(*difference, sep="\t", file=f)

    # Print the statistics
    logger.info("Printing the statistics")
    concordance, completion = print_statistics(
        sample_order, completion=completion, concordance=concordance, out=out,
    )

    # Print the concordance file
    logger.info("Printing the concordance file")
    print_concordance(concordance, out)

    return completion, concordance, sample_order


def get_acgt_geno_map(a1: str, a2: str) -> Dict[int, str]:
    """Decode the numerical genotypes {0, 1, 2} to ACGT genotypes.

    Note that for plink, 0 is homozugous a2, 1 is heterozygous and 2 is
    homozugous a1 in the BIM file.

    """
    return {
        0: f"{a2} {a2}",
        1: f"{a1} {a2}",
        2: f"{a1} {a1}",
        -1: "0 0",
    }


def read_duplicates(filename: str) -> Dict[str, Set[Tuple[str, str]]]:
    """Reads the duplicated samples from a file.

    Args:
        filename (str): the name of the file describing the duplicated samples.

    Returns:
        dict: a dictionary describing the duplicated samples.

    """
    duplicated_samples = collections.defaultdict(set)

    with open(filename) as f:
        for line in f:
            sample, fid, iid = line.rstrip().split()
            duplicated_samples[sample].add((fid, iid))

    return duplicated_samples


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exists with code 1.

    """
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such files")

    if not path.isfile(args.duplicated_samples):
        raise ProgramError(f"{args.duplicated_samples}: no such file")

    # Checking the concordance threshold
    if not 0 <= args.sample_concordance_threshold <= 1:
        raise ProgramError(
            f"sample-concordance-threshold: must be between 0 and 1 "
            f"(not {args.sample_concordance_threshold})",
        )

    # Checking the completion threshold
    if not 0 <= args.sample_completion_threshold <= 1:
        raise ProgramError(
            f"sample-completion-threshold: must be between 0 and 1 "
            f"(not {args.sample_completion_threshold})"
        )


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parses the command line options and arguments."""
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument(
        "-v", "--version", action="version",
        version=f"pyGenClean {SCRIPT_NAME} {__version__}",
    )

    # Adding the arguments and options
    add_args(parser)

    return parser.parse_args(argv)


def add_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments and options to the parser."""
    # The INPUT files
    group = parser.add_argument_group("Input File")
    group.add_argument(
        "--bfile", type=str, metavar="FILE", required=True,
        help="The input file prefix (will find the tped and tfam file by "
             "appending the prefix to .tped and .tfam, respectively.) The "
             "duplicated samples should have the same identification numbers "
             "(both family and individual ids.)",
    )
    group.add_argument(
        "--duplicated-samples", type=str, metavar="FILE", required=True,
        help="The file describing the duplicated samples. This file has no "
             "header and has two columns (sperated by a tabulation). The "
             "first column is the sample ID (e.g. NA12878). The second column "
             "is the IDs in the FAM file, separated by a comma (e.g. "
             "NA12878_1,NA12878_2).",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )
    group.add_argument(
        "--sample-completion-threshold", type=float, metavar="FLOAT",
        default=0.9,
        help="The completion threshold to consider a replicate when choosing "
             "the best replicates and for creating the composite samples. "
             "[default: %(default).1f]",
    )
    group.add_argument(
        "--sample-concordance-threshold", type=float, metavar="FLOAT",
        default=0.97,
        help="The concordance threshold to consider a replicate when choosing "
             "the best replicates and for creating the composite samples. "
             "[default: %(default).2f]",
    )

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default="dup_samples",
        help="The prefix of the output files. [default: %(default)s]",
    )
