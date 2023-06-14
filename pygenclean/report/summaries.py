"""Summaries for QC modules."""


import argparse
import re
from os import path
from pathlib import Path
from typing import Dict, Optional, Union

import jinja2
import numpy as np
import pandas as pd
from jinja2 import BaseLoader, Environment

from .utils import format_numbers


LABEL_RE = re.compile(r"[\/]")


class Summary():
    """Summary core object."""
    methods = ""
    results = ""

    """The summary core function."""
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args

    def get_results_template(self) -> jinja2.environment.Template:
        """Generate the Jinja2 template for the results."""
        return Environment(loader=BaseLoader).from_string(self.results)

    def get_methods_template(self) -> jinja2.environment.Template:
        """Generate the Jinja2 template for the methods."""
        return Environment(loader=BaseLoader).from_string(self.methods)

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        raise NotImplementedError()

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        raise NotImplementedError()

    def generate_results(self) -> str:
        """Generate the summary (results)."""
        template = self.get_results_template()
        return template.render(**self.get_results_information())

    def generate_methods(self) -> str:
        """Generate the methods summary."""
        template = self.get_methods_template()
        return template.render(**self.get_methods_information())


class SubsetSummary(Summary):
    """Subset summary."""
    methods = (
        "Subsets {{ subset_type }}{{ ' (' + reason + ')' if reason }} using "
        "_Plink_."
    )

    results = """
### Subset{{ " (" + reason + ")" if reason }}

{% if nb_markers %}
The file for marker exclusion contained {{ "{:,d}".format(nb_markers) }}
markers to {{ marker_excl_type }}. Out of a total of
{{ "{:,d}".format(nb_markers_before) }} markers,
{{ "{:,d}".format(nb_markers_after) }} remained
({{ "{:,d}".format(nb_markers_before - nb_markers_after) }} excluded).
{% endif %}

{% if nb_samples %}
The file for sample exclusion contained {{ "{:,d}".format(nb_samples) }}
samples to {{ sample_excl_type }}. Out of a total of
{{ "{:,d}".format(nb_samples_before) }} samples,
{{ "{:,d}".format(nb_samples_after) }} remained
({{ "{:,d}".format(nb_samples_before - nb_samples_after) }} excluded).
{% endif %}
"""

    def get_marker_info(self) -> Dict[str, Optional[int]]:
        """Get marker information."""
        nb_markers = None
        nb_markers_before = None
        nb_markers_after = None
        exclusion_type = None
        if self.args.extract or self.args.exclude:
            # Before
            with open(self.args.bfile + ".bim") as f:
                nb_markers_before = len(f.read().splitlines())

            # After
            with open(self.args.out + ".bim") as f:
                nb_markers_after = len(f.read().splitlines())

            # Exclusion
            exclusion_file = self.args.extract if self.args.extract \
                else self.args.exclude
            with open(exclusion_file) as f:
                nb_markers = len(f.read().splitlines())

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
            # Before
            with open(self.args.bfile + ".fam") as f:
                nb_samples_before = len(f.read().splitlines())

            # After
            with open(self.args.out + ".fam") as f:
                nb_samples_after = len(f.read().splitlines())

            # Exclusion
            exclusion_file = self.args.remove if self.args.remove \
                else self.args.keep
            with open(exclusion_file) as f:
                nb_samples = len(f.read().splitlines())

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
        return {
            "reason": self.args.reason,
            **self.get_marker_info(),
            **self.get_sample_info(),
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

    results = """
### Sex check

Using $F$ thresholds of {{ male_f }} and {{ female_f }} for males and females,
respectively, {{ "{:,d}".format(nb_problems) }}
sample{{ "s" if nb_problems > 1 }} had sex problems according to _Plink_.
{%- if nb_problems > 1 %}
@tbl-{{ label_prefix }}-results summarizes the sex problems encountered during
the analysis.
{%- endif -%}
{%- if figure_intensities %}
@fig-{{ label_prefix }}-intensities shows the $y$ intensities versus the $x$
intensities for each samples. Problematic samples are shown using triangles.
{%- endif -%}
{%- if figure_baf_lrr|length > 0 -%}
{%- if figure_baf_lrr|length == 1 %}
@fig-{{ label_prefix }}-baf_lrr-{{ figure_baf_lrr[0][1] }} shows
{%- else %}
@fig-{{ label_prefix }}-baf_lrr-{{ figure_baf_lrr[0][1] }} to
@fig-{{ label_prefix }}-baf_lrr-{{ figure_baf_lrr[-1][1] }} show
{%- endif %}
the log R ratio and the B allele frequency versus the position on chromosome X
and Y for the problematic samples.
{% endif %}

{% if nb_problems > 1 %}
{{ table }}

: Summarization of the gender problems encountered during Plink's analysis. HET
is the heterozygosity rate on the X chromosome. NOCALL is the percentage of no
calls on the Y chromosome. {{ "{#" }}tbl-{{ label_prefix }}-results}

{% if figure_intensities %}
![
    Gender check using Plink. Mean $x$ and $y$ intensities are shown for each
    sample. Males are shown in blue, and females in red. Triangles show
    problematic samples (green for males, mauve for females). Unknown gender
    are shown in gray.
]({{ figure_intensities }}){{ "{#" }}fig-{{ label_prefix }}-intensities}
{% endif %}

{% if figure_baf_lrr|length > 0 %}
{% for figure_path, sample_id in figure_baf_lrr %}
![
    Plots showing the log R ratio and the B allele frequency for chromosome X
    and Y (on the left and right, respectively) for sample {{ sample_id }}.
]({{ figure_path }}){{ "{#" }}fig-{{ label_prefix }}-baf_lrr-{{ sample_id }}}
{% endfor %}
{% endif %}
{% endif %}
"""

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

        return {
            "male_f": self.args.male_f,
            "female_f": self.args.female_f,
            "nb_problems": nb_problems,
            "table": sex_problems.to_markdown(index=False),
            "figure_intensities": figure_intensities,
            "figure_baf_lrr": list(zip(baf_lrr_figures, baf_lrr_samples)),
            "label_prefix": LABEL_RE.sub("-", self.args.out),
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

    results = """
### No calls and heterozygous only markers

After scrutiny, {{ "{:,d}".format(all_failed) }}
marker{{ "s" if all_failed > 1 }} {{ "was" if all_failed == 1 else "were" }}
excluded from the dataset because of a call rate of 0. Also,
{{ "{:,d}".format(all_hetero) }} marker{{ "s" if all_hetero > 1 }}
{{ "was" if all_hetero == 1 else "was" }} excluded from the dataset because all
samples were heterozygous (excluding the mitochondrial chromosome).
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # All failed markers
        with open(self.args.out + ".all_failed") as f:
            all_failed = len(f.read().splitlines())

        # All heterozygous markers
        with open(self.args.out + ".all_hetero") as f:
            all_hetero = len(f.read().splitlines())

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

    results = """
### Sample call rate

Using a `mind` threshold of {{ mind }} (_i.e._ keeping only samples with a
missing rate $\\leq{{ mind }}$), {{ nb_samples }}
sample{{ "s" if nb_samples > 1 }} {{ "was" if nb_samples == 1 else "were" }}
excluded from the dataset.
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        with open(self.args.bfile + ".fam") as f:
            nb_before = len(f.read().splitlines())

        with open(self.args.out + ".fam") as f:
            nb_after = len(f.read().splitlines())

        return {
            "mind": self.args.mind,
            "nb_samples": nb_before - nb_after,
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

    results = """
### Marker call rate

Using a `geno` threshold of {{ geno }} (_i.e._ keeping only markers with a
missing rate $\\leq{{ geno }}$), {{ nb_markers }}
marker{{ "s" if nb_markers > 1 }} {{ "was" if nb_markers == 1 else "were" }}
excluded from the dataset.
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        with open(self.args.bfile + ".bim") as f:
            nb_before = len(f.read().splitlines())

        with open(self.args.out + ".bim") as f:
            nb_after = len(f.read().splitlines())

        return {
            "geno": self.args.geno,
            "nb_markers": nb_before - nb_after,
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

    results = """
### Related samples

{% if not_enough -%}

There were not enough markers to perform the analysis
(_i.e._ {{ "{:,d}".format(nb_markers) }} markers is less than the required
{{ "{:,d}".format(nb_markers_required) }} markers).

{%- elif genome_only -%}

The _genome_ file was created using Plink and {{ "{:,d}".format(nb_markers) }}
markers.

{%- else -%}

According to _Plink_ relatedness analysis (using
{{ "{:,d}".format(nb_markers) }} markers),
{{ "{:,d}".format(nb_unique_samples) }} unique samples were related
to at least one other sample. A total of {{ "{:,d}".format(nb_discarded) }}
samples were randomly selected for downstream exclusion from the dataset.
{% if figure_z1 -%}
@fig-{{ label_prefix }}-z1 shows $Z_1$ versus $IBS2_{ratio}^\\ast$ for all
related samples found by _Plink_.
{% endif -%}
{% if figure_z2 -%}
@fig-{{ label_prefix }}-z2 shows $Z_2$ versus $IBS2_{ratio}^\\ast$ for all
related samples found by _Plink_.
{% endif %}

{% if figure_z1 %}
![
    $Z_1$ versus $IBS2_{ratio}^\\ast$ for all related samples found by _Plink_
    in the IBS analysis.
]({{ figure_z1 }}){{ "{#" }}fig-{{ label_prefix }}-z1}
{% endif %}

{% if figure_z2 %}
![
    $Z_2$ versus $IBS2_{ratio}^\\ast$ for all related samples found by _Plink_
    in the IBS analysis.
]({{ figure_z2 }}){{ "{#" }}fig-{{ label_prefix }}-z2}
{% endif %}

{%- endif %}
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Getting the number of pruned markers
        with open(self.args.out + ".pruned_data.bim") as f:
            nb_markers = len(f.read().splitlines())

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
        with open(self.args.out + ".discarded_related_individuals") as f:
            nb_discarded = len(f.read().splitlines())

        # Checking if we have z1 figure
        if path.isfile(self.args.out + ".related_individuals_z1.png"):
            figure_z1 = self.args.out + ".related_individuals_z1.png"

        # Checking if we have z2 figure
        if path.isfile(self.args.out + ".related_individuals_z2.png"):
            figure_z2 = self.args.out + ".related_individuals_z2.png"

        return {
            "nb_markers": nb_markers,
            "nb_unique_samples": len(related),
            "nb_discarded": nb_discarded,
            "figure_z1": figure_z1,
            "figure_z2": figure_z2,
            "label_prefix": LABEL_RE.sub("-", self.args.out),
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

    results = """
### Contamination

A total of {{ "{:,d}".format(nb_samples) }} sample{{ "s" if nb_samples > 1 }}
{{ "were" if nb_samples > 1 else "was" }} analyzed for contamination using
_bafRegress_. The analysis was performed using
{{ "{:,d}".format(nb_autosomal) }} autosomal
marker{{ "s" if nb_autosomal > 1 }}. Using a threshold of 0.01,
{{ "{:,d}".format(nb_contaminated) }} sample{{ "s" if nb_contaminated > 1 }}
{{ "were" if nb_contaminated > 1 else "was" }} estimated to be contaminated.
{%- if nb_contaminated > 1 %}
@tbl-{{ label_prefix }}-results lists all the samples that were estimated to be
contaminated (_i.e._ with an estimate $>0.01$)

{{ table }}

: List of all possible contaminated samples (_i.e._ with an estimate computed
by _bafRegress_ $>0.01$). {{ "{#" }}tbl-{{ label_prefix }}-results}
{% endif -%}
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Counting the number of autosomal markers
        nb_autosomal = 0
        with open(self.args.out + ".to_extract") as f:
            for _ in f:
                nb_autosomal += 1

        # Reading the output, and finding the contaminated samples
        df = pd.read_csv(self.args.out + ".bafRegress", sep="\t")
        contaminated = df.estimate > 0.01
        nb_contaminated = np.count_nonzero(contaminated)

        return {
            "nb_samples": df.shape[0],
            "nb_autosomal": nb_autosomal,
            "nb_contaminated": nb_contaminated,
            "table": df.loc[contaminated, :].to_markdown(index=False),
            "label_prefix": LABEL_RE.sub("-", self.args.out),
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

    results = """
### Plate bias

After performing the plate bias analysis using _Plink_, a total of
{{ "{:,d}".format(nb_significant) }} unique
marker{{ "s" if nb_significant != 1 }} had a significant result
(_i.e._ a value less than ${{ p_threshold }}$).
{%- if nb_significant > 0 %}
@tbl-{{ label_prefix }}-results summarizes the plate bias results.

{{ table }}

: Summary of the plate bias analysis performed by _Plink_. For each plate, the
number of significant markers is shown (threshold of ${{ p_threshold }}$). The
plates are sorted according to the total number of significant
results. {{ "{#" }}tbl-{{ label_prefix }}-results}
{% endif %}
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Reading the significant markers with the plate information
        df = pd.read_csv(self.args.out + ".significant_markers.summary.tsv",
                         sep="\t")
        return {
            "nb_significant": df.shape[0],
            "p_threshold": format_numbers(self.args.p_filter),
            "table": df.plate.value_counts()
                       .sort_values(ascending=False)
                       .to_markdown(),
            "label_prefix": LABEL_RE.sub("-", self.args.out),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}


class EthnicitySummary(Summary):
    """Plate bias summary."""
    methods = (
        "Checks sample's ethnicity using reference populations and IBS. The "
        "script uses pairwise IBS matrix as a distance metric to identify "
        "cryptic relatedness among samples and sample outliers by "
        "multidimensional scaling (MDS)."
    )

    results = """
### Ethnicity

{% if skip_ref_pops -%}

Principal components analysis was performed using
{{ "{:,d}".format(nb_markers) }} marker{{ "s" if nb_markers > 1 }} on the study
dataset only.

{%- else -%}

Using {{ "{:,d}".format(nb_markers) }} marker{{ "s" if nb_markers > 1 }} and a
multiplier of {{ multiplier }}, there was a total of
{{ "{:,d}".format(nb_outliers) }} outlier{{ "s" if nb_outliers > 1 }} of the
{{ outliers_of }} population.
{% if outlier_figure -%}
@fig-{{ label_prefix }}-outliers shows the first two principal components
of the MDS analysis, where outliers of the {{ outliers_of }} population are
shown in grey.
{% endif -%}
{% if scree_figure -%}
@fig-{{ label_prefix }}-scree shows the scree plot for the principal components
of the MDS analysis.
{% endif %}

{% if outlier_figure %}
![
    MDS plots showing the first two principal components of the source dataset
    with the reference panels. The outliers of the {{ outliers_of }} population
    are shown in grey, while samples of the source dataset that resemble the
    {{ outliers_of }} population are shown in orange. A multiplier of
    {{ multiplier }} was used to find the {{ "{:,d}".format(nb_outliers) }}
    outlier{{ "s" if nb_outliers > 1 }}.
]({{ outlier_figure }}){{ "{#" }}fig-{{ label_prefix }}-outliers}
{% endif %}

{% if scree_figure %}
![
    Scree plot for the principal components of the MDS analysis.
]({{ scree_figure }}){{ "{#" }}fig-{{ label_prefix }}-scree width='10cm'}
{% endif %}

{%- endif %}
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # The number of markers
        with open(self.args.out + ".ibs.pruned_data.bim") as f:
            nb_markers = len(f.read().splitlines())

        if self.args.skip_ref_pops:
            return {
                "nb_markers": nb_markers,
                "skip_ref_pops": self.args.skip_ref_pops,
            }

        # The number of outliers
        with open(self.args.out + ".outliers") as f:
            nb_outliers = len(f.read().splitlines())

        # The outlier figure
        outlier_figure = None
        if path.isfile(self.args.out + ".outliers.png"):
            outlier_figure = self.args.out + ".outliers.png"

        scree_figure = None
        if path.isfile(self.args.out + ".smartpca.scree_plot.png"):
            scree_figure = self.args.out + ".smartpca.scree_plot.png"

        return {
            "skip_ref_pops": self.args.skip_ref_pops,
            "nb_markers": nb_markers,
            "multiplier": self.args.multiplier,
            "outliers_of": self.args.outliers_of,
            "nb_outliers": nb_outliers,
            "outlier_figure": outlier_figure,
            "scree_figure": scree_figure,
            "label_prefix": LABEL_RE.sub("-", self.args.out),
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

    results = """
### Heterozygous haploid

After _Plink_'s heterozygous haploid analysis, a total of
{{ "{:,d}".format(nb_hetero_hap) }} genotype{{ "s" if nb_hetero_hap > 1}}
{{ "was" if nb_hetero_hap == 1 else "were" }} set to missing.
"""

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

    results = """
### Flag MAF

After computing minor allele frequencies (MAF) of all markers using _Plink_, a
total of {{ "{:,d}".format(nb_flagged) }} marker{{ "s" if nb_flagged > 1 }} had
a MAF of zero and were flagged (see file `flag_maf_0.list` for more
information).
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # the number of flagged markers
        with open(self.args.out + ".list") as f:
            nb_flagged = len(f.read().splitlines())

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

    results = """
### Flag Hardy-Weinberg

{% if flagged_markers|length == 0 %}

There were no markers to perform the analysis.

{% else %}

Markers which failed Hardy-Weinberg equilibrium test (using _Plink_) were
flagged.
{%- for p_threshold, nb_flagged, _ in flagged_markers %}
A total of {{ "{:,d}".format(nb_flagged) }} markers failed with a threshold of
${{ p_threshold }}$.
{%- endfor %}
For a total list of markers, check the files
{%- for _, _, filename in flagged_markers %}
{{ "and " if loop.last }}`{{ filename }}`{{ "," if loop.index < (loop.length - 1) }}
{%- endfor -%}, respectively.

{% endif %}
"""  # noqa: E501

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # Finding the files containing the flagged markers
        flag_files = Path(path.dirname(self.args.out))\
            .glob("flag_hw.snp_flag_threshold_*")

        p_threshold_re = re.compile(r"^flag_hw\.snp_flag_threshold_")

        flagged_markers = []
        for filename in flag_files:
            if "between" in filename.name:
                continue

            # The p threshold
            p_threshold = format_numbers(p_threshold_re.sub("", filename.name))

            # The number of flagged markers
            with open(filename) as f:
                nb_flagged = len(f.read().splitlines())

            flagged_markers.append((p_threshold, nb_flagged, filename.name))

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

    results = """
### Duplicated samples

A total of {{ "{:,d}".format(nb_dup_samples) }} duplicated
sample{{ "s" if nb_dup_samples > 1 }}
{{ "was" if nb_dup_samples == 1 else "were" }} found.
"""

    def get_results_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the results."""
        # The duplicated samples
        with open(self.args.duplicated_samples) as f:
            dup_samples = {line.split()[0] for line in f}

        return {
            "nb_dup_samples": len(dup_samples),
        }

    def get_methods_information(self) -> Dict[str, Optional[Union[str, int]]]:
        """Get the summary information for the methods."""
        return {}
