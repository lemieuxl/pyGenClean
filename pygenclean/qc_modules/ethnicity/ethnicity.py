"""Checks sample's ethnicity using reference populations"""


import argparse
import logging
from collections import defaultdict
from typing import Dict, KeysView, List, Optional, Set, Tuple

from ...error import ProgramError
from ...utils import flip_alleles
from ...utils import plink as plink_utils
from ...utils import timer
from ...utils.task import execute_external_command
from ...version import pygenclean_version as __version__
from ..related_samples import related_samples
from . import find_outliers, plot_eigenvalues, plot_mds


SCRIPT_NAME = "ethnicity"
DESCRIPTION = "Checks sample's ethnicity using reference populations."
DEFAULT_OUT = "ethnicity"


logger = logging.getLogger(__name__)


@timer(logger)
def main(args: Optional[argparse.Namespace] = None,
         argv: Optional[List[str]] = None) -> None:
    """Checks sample's ethnicity using reference populations.

    Args:
        args (argparse.Namespace): the arguments and options.
        argv (list): the argument as a list.

    These are the steps of this module:

    1.  Prints the options.
    2.  Finds the overlapping markers between the three reference panels and
        the source panel (:py:func:`findOverlappingSNPsWithReference`).
    3.  Extract the required markers from all the data sets
        (:py:func:`extractSNPs`).
    4.  Renames the reference panel's marker names to that they are the same as
        the source panel (for all populations) (:py:func:`renameSNPs`).
    5.  Combines the three reference panels together
        (:py:func:`combinePlinkBinaryFiles`).
    6.  Compute the frequency of all the markers from the reference and the
        source panels (:py:func:`computeFrequency`).
    7.  Finds the markers to flip in the reference panel (when compared to the
        source panel) (:py:func:`findFlippedSNPs`).
    8.  Excludes the unflippable markers from the reference and the source
        panels (:py:func:`excludeSNPs`).
    9.  Flips the markers that need flipping in their reference panel
        (:py:func:`flipSNPs`).
    10. Combines the reference and the source panels
        (:py:func:`combinePlinkBinaryFiles`).
    11. Runs part of :py:mod:`pyGenClean.RelatedSamples.find_related_samples`
        on the combined data set (:py:func:`runRelatedness`).
    12. Creates the ``mds`` file from the combined data set and the result of
        previous step (:py:func:`createMDSFile`).
    13. Creates the population file (:py:func:`createPopulationFile`).
    14. Plots the ``mds`` values (:py:func:`plotMDS`).
    15. Finds the outliers of a given reference population
        (:py:func:`find_the_outliers`).
    16. If required, computes the Eigenvalues using smartpca
        (:py:func:`compute_eigenvalues`).
    17. If required, creates a scree plot from smartpca resutls
        (:py:func:`create_scree_plot`).

    """
    if args is None:
        args = parse_args(argv)
    check_args(args)

    logger.info("%s", DESCRIPTION)

    ref_pop_names = ("CEU", "YRI", "JPT-CHB")

    bfile_for_genome = args.bfile
    if not args.skip_ref_pops:
        # The reference files information
        ref_pop_prefixes = (args.ceu_bfile, args.yri_bfile, args.jpt_chb_bfile)

        # Find overlap with the reference file
        logger.info("Finding overlapping SNPs between reference and "
                    "source panels")
        find_overlapping_snps_with_ref(
            prefix=args.bfile,
            ref_prefixes=ref_pop_prefixes,
            ref_pop_names=ref_pop_names,
            out_prefix=args.out,
        )

        # Extract the required SNPs using Plink (reference panels)
        logger.info("Extracting overlapping SNPs from the reference panels")
        for ref_name, ref_pop_prefix in zip(ref_pop_names, ref_pop_prefixes):
            logger.info("  - %s (%s)", ref_name, ref_pop_prefix)
            plink_utils.subset(
                bfile=ref_pop_prefix,
                out=args.out + f".reference_panel.{ref_name}",
                markers=args.out + f".{ref_name}_snp_to_extract",
                marker_subset_type="extract",
                use_original_plink=args.plink_107,
            )

        # Extract the required SNPs using Plink (source panel)
        logger.info("Extracting overlapping SNPs from the source panel")
        plink_utils.subset(
            bfile=args.bfile,
            out=args.out + ".source_panel",
            markers=args.out + ".source_snp_to_extract",
            marker_subset_type="extract",
            use_original_plink=args.plink_107,
        )

        # Renaming the reference file, so that the SNP names are the same
        for ref_name in ref_pop_names:
            plink_utils.rename_markers(
                bfile=args.out + f".reference_panel.{ref_name}",
                rename=args.out + f".{ref_name}_update_names",
            )

        # Combining the reference panel
        logger.info("Combining the reference panels")
        plink_utils.merge_files(
            prefixes=[
                args.out + f".reference_panel.{pop_name}"
                for pop_name in ref_pop_names
            ],
            out=args.out + ".reference_panel.ALL",
            use_original_plink=args.plink_107,
        )

        # Computing the frequency (reference panel)
        logger.info("Computing reference panel frequencies")
        plink_utils.compute_freq(
            bfile=args.out + ".reference_panel.ALL",
            out=args.out + ".reference_panel.ALL.frequency",
            use_original_plink=args.plink_107,
        )

        # Computing the frequency (source panel)
        logger.info("Computing source panel frequencies")
        plink_utils.compute_freq(
            bfile=args.out + ".source_panel",
            out=args.out + ".source_panel.frequency",
            use_original_plink=args.plink_107,
        )

        # Finding the SNPs to flip and flip them in the reference panel
        logger.info("Finding SNPs to flip or to exclude from reference panel")
        find_flipped_snps(
            frq_f1=args.out + ".reference_panel.ALL.frequency.frq",
            frq_f2=args.out + ".source_panel.frequency.frq",
            out=args.out,
        )

        # Excluding SNPs (reference panel)
        logger.info("Excluding SNPs from reference panel")
        plink_utils.subset(
            bfile=args.out + ".reference_panel.ALL",
            out=args.out + ".reference_panel.ALL.cleaned",
            markers=args.out + ".snp_to_remove",
            marker_subset_type="exclude",
            use_original_plink=args.plink_107,
        )

        # Excluding SNPs (source panel)
        logger.info("Excluding SNPs from source panel")
        plink_utils.subset(
            bfile=args.out + ".source_panel",
            out=args.out + ".source_panel.cleaned",
            markers=args.out + ".snp_to_remove",
            marker_subset_type="exclude",
            use_original_plink=args.plink_107,
        )

        # Flipping the SNP that need to be flip in the reference
        logger.info("Flipping SNPs in reference panel")
        plink_utils.flip_markers(
            prefix=args.out + ".reference_panel.ALL.cleaned",
            out=args.out + ".reference_panel.ALL.cleaned.flipped",
            markers=args.out + ".snp_to_flip_in_reference",
            use_original_plink=args.plink_107,
        )

        # Combining the reference panel
        logger.info("Combining reference and source panels")
        plink_utils.merge_files(
            prefixes=[args.out + ".reference_panel.ALL.cleaned.flipped",
                      args.out + ".source_panel.cleaned"],
            out=args.out + ".final_dataset_for_genome",
            use_original_plink=args.plink_107,
        )

        # The bfile for the genome file
        bfile_for_genome = args.out + ".final_dataset_for_genome"

    # Runing the relatedness step
    logger.info("Creating the genome file using Plink")
    related_samples.main(
        argv=[
            "--bfile", bfile_for_genome,
            "--genome-only",
            "--min-nb-snp", str(args.min_nb_snp),
            "--maf", args.maf,
            "--out", args.out + ".ibs",
            "--nb-threads", str(args.nb_threads),
            "--indep-pairwise",
        ] + args.indep_pairwise,
    )

    # Creating the MDS file
    logger.info("Creating the MDS file using Plink")
    create_mds(
        bfile=related_samples.get_prefix_for_genome(args.out + ".ibs"),
        out=args.out + ".mds",
        nb_components=args.nb_components,
        genome_filename=args.out + ".ibs.genome.genome.gz",
        use_original_plink=args.plink_107,
    )

    if not args.skip_ref_pops:
        # Creating the population files
        logger.info("Creating a population file")
        create_population_file(
            fam_files=[
                args.out + f".reference_panel.{pop_name}.fam"
                for pop_name in ref_pop_names
            ] + [args.out + ".source_panel.fam"],
            labels=list(ref_pop_names) + ["SOURCE"],
            out=args.out + ".population_file",
        )

        # Plot the MDS value
        logger.info("Creating the MDS plot")
        plot_mds.main(
            argv=[
                "--file", args.out + ".mds.mds",
                "--population-file", args.out + ".population_file",
                "--format", args.plot_mds_format,
                "--xaxis", args.plot_mds_xaxis,
                "--yaxis", args.plot_mds_yaxis,
                "--title", args.plot_mds_title,
                "--out", args.out + ".mds",
            ]
        )

        # Finding the outliers
        logger.info("Finding the outliers")
        find_outliers.main(
            argv=[
                "--mds", args.out + ".mds.mds",
                "--population-file", args.out + ".population_file",
                "--outliers-of", args.outliers_of,
                "--multiplier", str(args.multiplier),
                "--format", args.find_outliers_format,
                "--xaxis", args.find_outliers_xaxis,
                "--yaxis", args.find_outliers_yaxis,
                "--out", args.out,
            ],
        )

    # De we need to create a scree plot?
    if args.create_scree_plot:
        # Computing the eigenvalues using smartpca
        logger.info("Computing eigenvalues")
        compute_eigenvalues(
            prefix=args.out + ".ibs.pruned_data",
            out=args.out + ".smartpca",
            nb_components=args.nb_components,
        )

        logger.info("Creating scree plot")
        plot_eigenvalues.main(
            argv=[
                "--evec", args.out + ".smartpca.evec.txt",
                "--out", args.out + ".smartpca.scree_plot.png",
                "--title", args.plot_eigen_title,
            ],
        )

    return {
        "usable_bfile": args.bfile,
    }


def compute_eigenvalues(prefix: str, out: str, nb_components: int) -> None:
    """Computes the Eigenvalues using smartpca from Eigensoft.

    Args:
        prefix (str): the prefix of the input files.
        out (str): the prefix of the output files.

    Creates a "parameter file" used by smartpca and runs it.

    """
    # First, we create the parameter file
    with open(out + ".parameters", "w") as f:
        print(f"genotypename:    {prefix}.bed", file=f)
        print(f"snpname:         {prefix}.bim", file=f)
        print(f"indivname:       {prefix}.fam", file=f)
        print(f"numoutevec:      {nb_components}", file=f)
        print(f"evecoutname:     {out}.evec.txt", file=f)
        print(f"evaloutname:     {out}.eval.txt", file=f)
        print("numoutlieriter:  0", file=f)
        print("altnormstyle:    NO", file=f)

    # Executing smartpca
    command = [
        "smartpca",
        "-p", out + ".parameters",
    ]
    execute_external_command(command)


def create_population_file(fam_files: List[str], labels: List[str],
                           out: str) -> None:
    """Creates a population file.

    Args:
        fam_files (list): the list of FAM files.
        labels (list): the list of labels (corresponding to the FAM files).
        out (str): the name of the output file.

    The ``fam_files`` is in reality a list of ``fam`` files composed of
    samples. For each of those ``fam`` files, there is a label associated with
    it (representing the name of the population).

    The output file consists of one row per sample, with the following three
    columns: the family ID, the individual ID and the population of each
    sample.

    """
    with open(out, "w") as f:
        for fam_file, label in zip(fam_files, labels):
            for row in plink_utils.parse_fam(fam_file):
                print(row.fid, row.iid, label, sep="\t", file=f)


def create_mds(bfile: str, out: str, genome_filename: str, nb_components: int,
               use_original_plink: bool = False) -> None:
    """Creates a MDS file using Plink.

    Args:
        bfile (str): the prefix of the input file.
        out (str): the prefix of the output file.
        genome_filename (str): the name of the ``genome`` file.
        nb_components (int): the number of component.
        use_original_plink (bool): plink1.9 or plink1.07.

    Using Plink, computes the MDS values for each individual using the
    ``inPrefix``, ``genomeFileName`` and the number of components. The results
    are save using the ``outPrefix`` prefix.

    """
    command = [
        "plink" if use_original_plink else "plink1.9",
        "--noweb",
        "--bfile", bfile,
        "--read-genome", genome_filename,
        "--cluster",
        "--mds-plot", str(nb_components),
        "--out", out,
    ]
    execute_external_command(command)


def find_flipped_snps(frq_f1: str, frq_f2: str, out: str) -> None:
    """Find flipped SNPs and flip them in the data.

    Args:
        frq_f1 (str): the name of the first frequency file
        frq_f2 (str): the name of the second frequency file.
        out (str): the prefix of the output files.

    By reading two frequency files (``frqFile1`` and ``frqFile2``), it finds a
    list of markers that need to be flipped so that the first file becomes
    comparable with the second one. Also finds marker that need to be removed.

    A marker needs to be flipped in one of the two data set if the two markers
    are not comparable (same minor allele), but become comparable if we flip
    one of them.

    A marker will be removed if it is all homozygous in at least one data set.
    It will also be removed if it's impossible to determine the phase of the
    marker (*e.g.* if the two alleles are ``A`` and ``T`` or ``C`` and ``G``).

    """
    alleles_f1: Dict[str, Set[str]] = {}
    alleles_f2: Dict[str, Set[str]] = {}
    for alleles_f, filename in zip((alleles_f1, alleles_f2), (frq_f1, frq_f2)):
        with open(filename, "r") as f:
            header = None
            for line in f:
                row = plink_utils.split_line(line)

                if header is None:
                    header = {name: i for i, name in enumerate(row)}

                    # Checking the columns
                    for column in ["SNP", "A1", "A2"]:
                        if column not in header:
                            raise ProgramError(
                                f"{filename}: no column named {column}",
                            )

                    continue

                # The alleles
                alleles = set([row[header["A1"]], row[header["A2"]]])

                # Removing "unknown" allele
                if "0" in alleles:
                    alleles.remove("0")

                alleles_f[row[header["SNP"]]] = alleles

    # Finding the SNPs to flip
    to_flip_f = open(out + ".snp_to_flip_in_reference", "w")

    to_remove_f = open(out + ".snp_to_remove", "w")

    for snp_name, alleles_1 in alleles_f1.items():
        alleles_2 = alleles_f2[snp_name]

        if (len(alleles_1) == 2) and (len(alleles_2) == 2):
            # Both are heterozygous
            if ({"A", "T"} == alleles_1) or ({"C", "G"} == alleles_1) or \
                    ({"A", "T"} == alleles_2) or ({"C", "G"} == alleles_2):
                # We can't flip those..., so we remove them
                print(snp_name, file=to_remove_f)

            else:
                if alleles_1 != alleles_2:
                    # Let's try the flip one
                    if flip_alleles(alleles_1) == alleles_2:
                        # We need to flip it
                        print(snp_name, file=to_flip_f)

                    else:
                        # Those SNP are discordant...
                        print(snp_name, file=to_remove_f)
        else:
            # We want to remove this SNP, because there is at least one
            # homozygous individual
            print(snp_name, file=to_remove_f)

    # Closing output files
    to_flip_f.close()
    to_remove_f.close()


def find_overlapping_snps_with_ref(
    prefix: str,
    ref_prefixes: Tuple[str, str, str],
    ref_pop_names: Tuple[str, str, str],
    out_prefix: str,
) -> None:
    """Find the overlapping SNPs in 4 different data sets.

    Args:
        prefix (str): the prefix source file.
        ref_prefixes (tuple): the prefix of the reference populations.
        ref_pop_names (tuple): the name of the reference populations.
        out_prefix (str): the prefix of the output files.

    It starts by reading the ``bim`` file of the source data set
    (``prefix.bim``). It finds all the markers (excluding the duplicated ones).
    Then it reads all of the reference population ``bim`` files
    (``referencePrefixes.bim``) and find all the markers that were found in the
    source data set.

    It creates three output files:

    * ``outPrefix.ref_snp_to_extract``: the name of the markers that needs to
      be extracted from the three reference panels.
    * ``outPrefix.source_snp_to_extract``: the name of the markers that needs
      to be extracted from the source panel.
    * ``outPrefix.update_names``: a file (readable by Plink) that will help in
      changing the names of the selected markers in the reference panels, so
      that they become comparable with the source panel.

    """
    # Reading the source BIM file
    source_snps = get_unique_markers(prefix)
    logger.info(
        "  - found %d unique markers in %s", len(source_snps), prefix + ".bim",
    )

    # Intersecting with the reference file
    ref_snps: Dict[str, Dict[Tuple[int, int], str]] = defaultdict(dict)
    for ref_name, ref_prefix in zip(ref_pop_names, ref_prefixes):
        ref_snps[ref_name] = get_unique_markers(
            ref_prefix, source_snps.keys(),
        )

        logger.info(
            "  - %d overlaps with %s", len(ref_snps[ref_name]), ref_prefix,
        )

    # Creating the intersect of the reference SNP
    ref_snp_intersection = set(ref_snps[ref_pop_names[0]].keys())
    for ref_name in ref_pop_names[1:]:
        ref_snp_intersection &= ref_snps[ref_name].keys()

    logger.info(
        "  - %d in common between reference panels", len(ref_snp_intersection),
    )

    # Printing the extraction and rename files
    with open(out_prefix + ".source_snp_to_extract", "w") as f, \
         open(out_prefix + ".CEU_snp_to_extract", "w") as list_ceu_f, \
         open(out_prefix + ".YRI_snp_to_extract", "w") as list_yri_f, \
         open(out_prefix + ".JPT-CHB_snp_to_extract", "w") as list_jpt_chb_f, \
         open(out_prefix + ".CEU_update_names", "w") as name_ceu_f, \
         open(out_prefix + ".YRI_update_names", "w") as name_yri_f, \
         open(out_prefix + ".JPT-CHB_update_names", "w") as name_jpt_chb_f:
        for name in ref_snp_intersection:
            # The extraction files
            print(source_snps[name], file=f)
            print(ref_snps["CEU"][name], file=list_ceu_f)
            print(ref_snps["YRI"][name], file=list_yri_f)
            print(ref_snps["JPT-CHB"][name], file=list_jpt_chb_f)

            # The rename files
            print(
                ref_snps["CEU"][name], source_snps[name],
                sep="\t", file=name_ceu_f,
            )
            print(
                ref_snps["YRI"][name], source_snps[name],
                sep="\t", file=name_yri_f,
            )
            print(
                ref_snps["JPT-CHB"][name], source_snps[name],
                sep="\t", file=name_jpt_chb_f,
            )


def get_unique_markers(
    prefix: str,
    subset: Optional[KeysView[Tuple[int, int]]] = None,
) -> Dict[Tuple[int, int], str]:
    """Gets the unique set of markers in a BIM file."""
    snps = {}
    duplicates = set()

    for row in plink_utils.parse_bim(prefix + ".bim"):
        marker_key = (row.chrom, row.pos)
        name = row.name

        if subset:
            if marker_key not in subset:
                continue

        # First time we saw it
        if marker_key not in snps:
            snps[marker_key] = name

        # It's a duplicate
        else:
            duplicates.add(marker_key)

    # Removing duplicates from the list
    for snp in duplicates:
        del snps[snp]

    return snps


def check_args(args: argparse.Namespace) -> None:
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :py:class:`sys.stderr` and the program exists with code 1.

    """
    # Check if we have the required files
    if not plink_utils.check_files(args.bfile):
        raise ProgramError(f"{args.bfile}: no such binary files")

    # Checking the reference population, if required
    if not args.skip_ref_pops:
        populations = ("ceu", "yri", "jpt-chb")
        prefixes = (args.ceu_bfile, args.yri_bfile, args.jpt_chb_bfile)
        for pop, prefix in zip(populations, prefixes):
            # Checking the option is present
            if prefix is None:
                raise ProgramError(
                    f"Missing {pop.upper()} population (--{pop}-bfile)",
                )

            # Checking the binary files are present
            if not plink_utils.check_files(prefix):
                raise ProgramError(f"{prefix}: no such binary files")

    # Check the indep-pairwise option
    # The two first must be int, the last one float
    try:
        _ = int(args.indep_pairwise[0])
        _ = int(args.indep_pairwise[1])
        _ = float(args.indep_pairwise[2])
    except ValueError as error:
        raise ProgramError("indep-pairwise: need INT INT FLOAT") from error

    # Check the maf value
    tmp_maf = None
    try:
        tmp_maf = float(args.maf)
    except ValueError as error:
        raise ProgramError(f"maf: must be a float, not {args.maf}") from error
    if (tmp_maf > 0.5) or (tmp_maf < 0.0):
        raise ProgramError("maf: must be between 0.0 and 0.5, not {args.maf}")

    # Check the number of component to compute
    if args.nb_components < 2:
        raise ProgramError(
            f"{args.nb_components} components: must be higher than 2",
        )

    # Check the minimum number of SNPs
    if args.min_nb_snp < 1:
        raise ProgramError(f"minimum {args.min_nb_snp} snp: must be above 1")


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
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively.",
    )
    group.add_argument(
        "--skip-ref-pops", action="store_true",
        help="Perform the MDS computation, but skip the three reference "
             "panels.",
    )
    group.add_argument(
        "--ceu-bfile", type=str, metavar="FILE",
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively) for the CEU population",
    )
    group.add_argument(
        "--yri-bfile", type=str, metavar="FILE",
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively) for the YRI population",
    )
    group.add_argument(
        "--jpt-chb-bfile", type=str, metavar="FILE",
        help="The input file prefix (will find the plink binary files by "
             "appending the prefix to the .bim, .bed and .fam files, "
             "respectively) for the JPT-CHB population",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--plink-1.07", dest="plink_107", action="store_true",
        help="Use original Plink (version 1.07)",
    )
    group.add_argument(
        "--nb-threads", type=int, metavar="N", default=1,
        help="The number of threads for this analysis (no effect when using "
             "plink 1.07). [%(default)d]",
    )
    group.add_argument(
        "--min-nb-snp", type=int, metavar="INT", default=8000,
        help="The minimum number of markers needed to compute IBS values. "
             "[Default: %(default)d]",
    )
    group.add_argument(
        "--indep-pairwise", type=str, metavar="STR", nargs=3,
        default=["50", "5", "0.1"],
        help="Three numbers: window size, window shift and the r2 threshold. "
             "[default: %(default)s]",
    )
    group.add_argument(
        "--maf", type=str, metavar="FLOAT", default="0.05",
        help="Restrict to SNPs with MAF >= threshold. [default: %(default)s]",
    )
    group.add_argument(
        "--nb-components", type=int, metavar="INT", default=10,
        help="The number of component to compute. [default: %(default)d]",
    )

    # The outlier options
    group = parser.add_argument_group("Outlier Options")
    find_outliers.add_graphical_options(group, prefix="find-outliers-")

    # The MDS plotting options
    group = parser.add_argument_group("MDS Plot Options")
    plot_mds.add_graphical_options(group, prefix="plot-mds-")

    # The Scree Plot options
    group = parser.add_argument_group("Scree Plot Options")
    group.add_argument(
        "--create-scree-plot", action="store_true",
        help="Computes Eigenvalues and creates a scree plot.",
    )
    plot_eigenvalues.add_graphical_options(group, prefix="plot-eigen-")

    # The OUTPUT files
    group = parser.add_argument_group("Output File")
    group.add_argument(
        "--out", type=str, metavar="FILE", default=DEFAULT_OUT,
        help="The prefix of the output files. [default: %(default)s]",
    )
