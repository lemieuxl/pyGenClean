"""Summaries for QC modules."""


import argparse
import re
from os import path
from pathlib import Path
from typing import Dict, Optional, Union

import jinja2
import numpy as np
import pandas as pd
from jinja2 import BaseLoader, Environment, PackageLoader

from ..utils import count_lines
from ..utils.plink import compare_bim, compare_fam
from .utils import format_numbers


LABEL_RE = re.compile(r"[\/]")


_TEMPLATES = Environment(loader=PackageLoader(__name__))


class Summary():
    """Summary core object."""
    methods = ""
    result_section_name = ""

    """The summary core function."""
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.summary_table_info = None
        self.label_prefix = LABEL_RE.sub("-", self.args.out)

    def get_results_template(self) -> jinja2.environment.Template:
        """Generate the Jinja2 template for the results."""
        return _TEMPLATES.get_template(self.__class__.__name__ + ".md")

    def get_methods_template(self) -> jinja2.environment.Template:
        """Generate the Jinja2 template for the methods."""
        return Environment(loader=BaseLoader).from_string(self.methods)

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        raise NotImplementedError()

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        raise NotImplementedError()

    def generate_results(self, **kwargs) -> str:
        """Generate the summary (results)."""
        return self.get_results_template().render(
            section_level=3,
            section_key="results",
            section_name=self.result_section_name,
            label_prefix=self.label_prefix,
            **kwargs,
            **self.get_results_information(),
        )

    def generate_methods(self) -> str:
        """Generate the methods summary."""
        template = self.get_methods_template()
        return template.render(**self.get_methods_information())

    def write_results(self, **kwargs) -> str:
        """Write the results."""
        filename = self.args.out + ".summary.qmd"
        with open(self.args.out + ".summary.qmd", "w") as f:
            print(self.generate_results(**kwargs), file=f)
        return filename


class SubsetSummary(Summary):
    """Subset summary."""
    methods = (
        "Subsets {{ subset_type }}{{ ' (' + reason + ')' if reason }} using "
        "_Plink_."
    )

    result_section_name = "Subset"

    def get_marker_info(self) -> Dict[str, Optional[int]]:
        """Get marker information."""
        nb_markers = None
        nb_markers_before = None
        nb_markers_after = None
        exclusion_type = None
        if self.args.extract or self.args.exclude:
            before, both, after = compare_bim(
                self.args.bfile + ".bim", self.args.out + ".bim",
            )
            assert len(after) == 0

            nb_markers_before = len(before) + len(both)
            nb_markers_after = len(both)

            # Exclusion
            exclusion_file = self.args.extract if self.args.extract \
                else self.args.exclude
            nb_markers = count_lines(exclusion_file)

            # Exclusion type
            exclusion_type = "extract" if self.args.extract else "exclude"

        return {
            "nb_markers": nb_markers,
            "nb_markers_before": nb_markers_before,
            "nb_markers_after": nb_markers_after,
            "marker_excl_type": exclusion_type,
        }

    def get_sample_info(self) -> Dict[str, Optional[int]]:
        """Get sample information."""
        nb_samples = None
        nb_samples_before = None
        nb_samples_after = None
        exclusion_type = None
        if self.args.keep or self.args.remove:
            before, both, after = compare_fam(
                self.args.bfile + ".fam", self.args.out + ".fam",
            )
            assert len(after) == 0

            nb_samples_before = len(before) + len(both)
            nb_samples_after = len(both)

            # Exclusion
            exclusion_file = self.args.remove if self.args.remove \
                else self.args.keep
            nb_samples = count_lines(exclusion_file)

            # Exclusion type
            exclusion_type = "keep" if self.args.keep else "remove"

        return {
            "nb_samples": nb_samples,
            "nb_samples_before": nb_samples_before,
            "nb_samples_after": nb_samples_after,
            "sample_excl_type": exclusion_type,
        }

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        summary_table_info = []

        # The marker information (if any)
        marker_info = self.get_marker_info()
        if marker_info["marker_excl_type"]:
            marker_label = (
                "excluded" if marker_info["marker_excl_type"] == "exclude"
                else "extracted"
            )
            marker_label = "Markers " + marker_label
            nb_markers = (
                marker_info["nb_markers_before"]
                - marker_info["nb_markers_after"]
            )

            summary_table_info.append((marker_label, nb_markers))

        # The sample information
        sample_info = self.get_sample_info()
        if sample_info["sample_excl_type"]:
            sample_label = (
                "kept" if sample_info["sample_excl_type"] == "keep"
                else "removed"
            )
            sample_label = "Samples " + sample_label
            nb_samples = (
                sample_info["nb_samples_before"]
                - sample_info["nb_samples_after"]
            )

            summary_table_info.append((sample_label, nb_samples))

        self.summary_table_info = tuple(summary_table_info)

        return {
            "reason": self.args.reason,
            **marker_info,
            **sample_info,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        subset_type = []
        if self.args.keep or self.args.remove:
            subset_type.append("samples")
        if self.args.exclude or self.args.extract:
            subset_type.append("markers")

        return {
            "subset_type": " and ".join(subset_type),
            "reason": self.args.reason,
        }


class SexCheckSummary(Summary):
    """Sexcheck summary."""
    methods = (
        "Check sample's genetic sex using _Plink_. The script identifies any "
        "individual with discrepancies between phenotype and genotype data "
        "for sex. Individuals with sex error are to be investigated."
    )

    result_section_name = "Sex check"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Reading samples with sex problems
        sex_problems = pd.read_csv(self.args.out + ".list_problem_sex",
                                   sep="\t")\
            .set_index(["FID", "IID"], verify_integrity=True)
        nb_problems = sex_problems.shape[0]

        # Adding the heterozygosity on chromosome 23
        heterozygosity = pd.read_csv(self.args.out + ".chr23.hetero.tsv",
                                     sep="\t")\
            .set_index(["FID", "IID"], verify_integrity=True)

        # Adding the no call frequency on chromosome 24
        no_call = pd.read_csv(self.args.out + ".chr24.no_call.tsv", sep="\t")\
            .set_index(["FID", "IID"], verify_integrity=True)

        # Merging
        sex_problems = sex_problems.assign(
            HET=heterozygosity.HETERO,
            NOCALL=no_call.F_MISS,
        ).reset_index()

        # Checking if we have a figure for sumarized intensities
        figure_intensities = None
        if path.isfile(self.args.out + ".png"):
            figure_intensities = self.args.out + ".png"

        # Checking if we have BAF and LRR figures
        baf_lrr_figures = list(Path(self.args.out + ".BAF_LRR").glob("*.png"))
        baf_lrr_samples = [
            re.match(r"sample_(\S+)_lrr_baf\.png$", file_path.name).group(1)
            for file_path in baf_lrr_figures
        ]

        no_genetic_sex = sex_problems.SNPSEX == 0
        self.summary_table_info = (
            ("No genetic sex", np.count_nonzero(no_genetic_sex)),
            ("Discordant sex", np.count_nonzero(~no_genetic_sex)),
        )

        return {
            "male_f": self.args.male_f,
            "female_f": self.args.female_f,
            "nb_problems": nb_problems,
            "table": sex_problems.to_markdown(index=False, floatfmt=".4f"),
            "figure_intensities": figure_intensities,
            "figure_baf_lrr": list(zip(baf_lrr_figures, baf_lrr_samples)),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class NoCallHeteroSummary(Summary):
    """No call and heterozygotes summary."""
    methods = (
        "Removes completely failed markers (no calls) or markers with only "
        "heterozygous genotypes (excluding the mitochondrial chromosome)."
    )

    result_section_name = "No calls and heterozygous only markers"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # All failed markers
        all_failed = count_lines(self.args.out + ".all_failed")

        # All heterozygous markers
        all_hetero = count_lines(self.args.out + ".all_hetero")

        self.summary_table_info = (
            ("All failed", all_failed),
            ("All heterozygous", all_hetero),
        )

        return {
            "all_failed": all_failed,
            "all_hetero": all_hetero,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class SampleCallRateSummary(Summary):
    """Sample call rate summary."""
    methods = (
        "Computes sample call rate using Plink (`mind` = {{ mind }}). The "
        "script identifies poorly performing DNA samples with a genome-wide "
        "genotyping success rate $<{{ (1 - mind) * 100 }}\\%$."
    )

    result_section_name = "Sample call rate"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        before, _, after = compare_fam(
            self.args.bfile + ".fam", self.args.out + ".fam",
        )
        assert len(after) == 0

        self.summary_table_info = (
            ("mind=" + str(self.args.mind), len(before)),
        )

        return {
            "mind": self.args.mind,
            "nb_samples": len(before),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {
            "mind": self.args.mind,
        }


class MarkerCallRateSummary(Summary):
    """Marker call rate summary."""
    methods = (
        "Computes marker call rate using Plink (`geno` = {{ geno }}). The "
        "script identifies poorly performing markers genotyping success rate "
        "$<{{ (1 - geno) * 100 }}\\%$."
    )

    result_section_name = "Marker call rate"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        before, _, after = compare_bim(
            self.args.bfile + ".bim", self.args.out + ".bim"
        )
        assert len(after) == 0

        self.summary_table_info = (
            ("geno=" + str(self.args.geno), len(before)),
        )

        return {
            "geno": self.args.geno,
            "nb_markers": len(before),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {
            "geno": self.args.geno,
        }


class RelatedSamplesSummary(Summary):
    """Related samples summary."""
    methods = (
        "Finds related samples according to IBS values. The script conducts "
        "close familial relationship checks with pairwise IBD. It flags and "
        "removes all but one pair member of samples duplicates "
        "($IBS2^{\\ast}_{ratio} > {{ ibs2_ratio }}$) based on a selection of "
        "uncorrelated markers ($r^2 < {{ r2_threshold }}$)."
    )

    result_section_name = "Related samples"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Getting the number of pruned markers
        nb_markers = count_lines(self.args.out + ".pruned_data.bim")

        # Did we have enough?
        if nb_markers < self.args.min_nb_snp:
            return {
                "nb_markers": nb_markers,
                "nb_markers_required": self.args.min_nb_snp,
                "not_enough": True,
            }

        if self.args.genome_only:
            return {
                "genome_only": True,
                "nb_markers": nb_markers,
            }

        # Reading the related sample file
        related_samples = pd.read_csv(self.args.out + ".related_individuals",
                                      sep="\t")

        # Counting the number of related samples
        related = set()
        for _, row in related_samples.iterrows():
            related.add((row.FID1, row.IID1))
            related.add((row.FID2, row.IID2))

        # Reading the number of discated samples
        nb_discarded = count_lines(
            self.args.out + ".discarded_related_individuals",
        )

        # Checking if we have z1 figure
        if path.isfile(self.args.out + ".related_individuals_z1.png"):
            figure_z1 = self.args.out + ".related_individuals_z1.png"

        # Checking if we have z2 figure
        if path.isfile(self.args.out + ".related_individuals_z2.png"):
            figure_z2 = self.args.out + ".related_individuals_z2.png"

        self.summary_table_info = (
            ("Markers used for IBS", nb_markers),
            ("Unique related samples", len(related)),
        )

        # Reading the merged related samples
        merged_related_samples = pd.read_csv(
            self.args.out + ".merged_related_individuals",
            sep="\t",
            usecols=["index", "FID1", "IID1", "FID2", "IID2", "status"],
        )
        table = None
        if merged_related_samples.shape[0]:
            table = merged_related_samples.to_markdown(index=False)

        return {
            "nb_markers": nb_markers,
            "nb_unique_samples": len(related),
            "nb_discarded": nb_discarded,
            "figure_z1": figure_z1,
            "figure_z2": figure_z2,
            "table": table,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {
            "ibs2_ratio": self.args.ibs2_ratio,
            "r2_threshold": self.args.indep_pairwise[-1],
        }


class ContaminationSummary(Summary):
    """Contamination summary."""
    methods = (
        "Check _BAF_ and _LogR_ ratio for data contamination. The script "
        "search for sample contamination using the `bafRegress.py` software."
    )

    result_section_name = "Contamination"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Counting the number of autosomal markers
        nb_autosomal = 0
        with open(self.args.out + ".to_extract") as f:
            for _ in f:
                nb_autosomal += 1

        # Reading the output, and finding the contaminated samples
        df = pd.read_csv(self.args.out + ".bafRegress", sep="\t")
        contaminated = df.estimate > self.args.estimate_threshold
        nb_contaminated = np.count_nonzero(contaminated)

        self.summary_table_info = (
            ("Possibly contaminated", nb_contaminated),
        )

        return {
            "nb_samples": df.shape[0],
            "nb_autosomal": nb_autosomal,
            "nb_contaminated": nb_contaminated,
            "threshold": self.args.estimate_threshold,
            "table": df.loc[contaminated, :].to_markdown(
                index=False, floatfmt=".4f",
            ),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class PlateBiasSummary(Summary):
    """Plate bias summary."""
    methods = (
        "Checks for plate bias. The script identifies markers with plate "
        "bias, based on the plates used to dilute DNA samples."
    )

    result_section_name = "Plate bias"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Reading the significant markers with the plate information
        df = pd.read_csv(self.args.out + ".significant_markers.summary.tsv",
                         sep="\t")

        p_threshold = format_numbers(self.args.p_filter)

        self.summary_table_info = (
            (f"Marker with plate bias ($p<{p_threshold}$)", df.shape[0]),
        )

        return {
            "nb_significant": df.shape[0],
            "p_threshold": p_threshold,
            "table": df.plate.value_counts()
                       .sort_values(ascending=False)
                       .to_markdown(),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class AncestrySummary(Summary):
    """Plate bias summary."""
    methods = (
        "Checks sample's ancestry using reference populations and IBS. The "
        "script uses pairwise IBS matrix as a distance metric to identify "
        "cryptic relatedness among samples and sample outliers by "
        "multidimensional scaling (MDS)."
    )

    result_section_name = "Ancestry"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # The number of markers
        nb_markers = count_lines(self.args.out + ".ibs.pruned_data.bim")

        # Screen plot might be there if skipping reference populations
        scree_figure = None
        if path.isfile(self.args.out + ".smartpca.scree_plot.png"):
            scree_figure = self.args.out + ".smartpca.scree_plot.png"

        if self.args.skip_ref_pops:
            self.summary_table_info = (
                ("Markers used for MDS", nb_markers),
            )

            return {
                "nb_markers": nb_markers,
                "skip_ref_pops": self.args.skip_ref_pops,
                "scree_figure": scree_figure,
            }

        # The number of outliers
        nb_outliers = count_lines(self.args.out + ".outliers")

        # The outlier figure
        outlier_figure = None
        if path.isfile(self.args.out + ".outliers.png"):
            outlier_figure = self.args.out + ".outliers.png"

        self.summary_table_info = (
            ("Markers used for MDS", nb_markers),
            (f"{self.args.outliers_of} outliers", nb_outliers),
        )

        return {
            "skip_ref_pops": self.args.skip_ref_pops,
            "nb_markers": nb_markers,
            "multiplier": self.args.multiplier,
            "outliers_of": self.args.outliers_of,
            "nb_outliers": nb_outliers,
            "outlier_figure": outlier_figure,
            "scree_figure": scree_figure,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class HeteroHapSummary(Summary):
    """Heterozygous haploid summary."""
    methods = (
        "Removes heterozygous haploid genotypes. The script sets heterozygous "
        "haploid genotypes to missing."
    )

    result_section_name = "Heterozygous haploid"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Regular expression to gather the number of heterozygous haploid
        plink_re = re.compile(r"Warning: (\d+) het. haploid genotypes present")
        if self.args.plink_107:
            plink_re = re.compile(
                r"(\d+) heterozygous haploid genotypes; set to missing",
            )

        with open(self.args.out + ".log") as f:
            nb_hetero_hap = int(plink_re.search(f.read()).group(1))

        self.summary_table_info = (
            ("Heterozygous haploid genotypes", nb_hetero_hap),
        )

        return {
            "nb_hetero_hap": nb_hetero_hap,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class FlagMafSummary(Summary):
    """Flag MAF summary."""
    methods = (
        "Flags markers with MAF of 0. The script flags markers with a minor "
        "allele frequency of 0 (_i.e._ monomorphic markers)."
    )

    result_section_name = "Flag MAF"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # the number of flagged markers
        nb_flagged = count_lines(self.args.out + ".list")

        self.summary_table_info = (
            ("Markers with a MAF of 0", nb_flagged),
        )

        return {
            "nb_flagged": nb_flagged,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class FlagHwSummary(Summary):
    """Flag HW summary."""
    methods = (
        "Flags markers with Hardy-Weinberg disequilibrium. The script tests "
        "for Hardy-Weinberg equilibrium for each marker (using an exact "
        "test). It adjusts for multiple testing using Bonferroni."
    )

    result_section_name = "Flag Hardy-Weinberg"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Finding the files containing the flagged markers
        flag_files = Path(path.dirname(self.args.out))\
            .glob("flag_hw.snp_flag_threshold_*")

        p_threshold_re = re.compile(r"^flag_hw\.snp_flag_threshold_")

        summary_table_info = []

        flagged_markers = []
        for filename in flag_files:
            if "between" in filename.name:
                continue

            # The p threshold
            p_threshold = float(p_threshold_re.sub("", filename.name))
            p_threshold = format_numbers(f"{p_threshold:.4g}")

            # The number of flagged markers
            nb_flagged = count_lines(filename)

            flagged_markers.append((p_threshold, nb_flagged, filename.name))

            summary_table_info.append((f"$p<{p_threshold}$", nb_flagged))

        self.summary_table_info = tuple(summary_table_info)

        return {
            "flagged_markers": flagged_markers,
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class DuplicatedSamplesSummary(Summary):
    """Duplicated samples summary."""
    methods = (
        "Extracts and merges duplicated samples. The script evaluates "
        "concordance and completion rate. If the thresholds are met, the "
        "script merges and completes the genotypes."
    )

    result_section_name = "Duplicated samples"

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # The duplicated samples
        with open(self.args.duplicated_samples) as f:
            dup_samples = {line.split()[0] for line in f}

        self.summary_table_info = (
            ("Number of duplicated samples", len(dup_samples)),
        )

        return {
            "nb_dup_samples": len(dup_samples),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}
