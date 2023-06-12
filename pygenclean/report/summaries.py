"""Summaries for QC modules"""


import argparse
import re
from os import path
from pathlib import Path

import pandas as pd
from jinja2 import BaseLoader, Environment


class SexCheckSummary():
    """Sexcheck summary."""
    def __init__(self, args: argparse.Namespace):
        self.args = args

        # The summary
        summary = """
### Sex check

Using $F$ thresholds of {{ male_f }} and {{ female_f }} for males and females,
respectively, {{ "{:,d}".format(nb_problems) }}
sample{{ "s" if nb_problems > 1 }} had sex problems according to _Plink_.
{%- if nb_problems > 1 %}
@tbl-sexcheck-results summarizes the sex problems encountered during the
analysis.
{%- endif -%}
{%- if figure_intensities %}
@fig-sexcheck-intensities shows the $y$ intensities versus the $x$ intensities
for each samples. Problematic samples are shown using triangles.
{%- endif -%}
{%- if figure_baf_lrr|length > 0 -%}
{%- if figure_baf_lrr|length == 1 %}
@fig-baf_lrr-{{ figure_baf_lrr[0][1] }} shows
{%- else %}
@fig-baf_lrr-{{ figure_baf_lrr[0][1] }} to
@fig-baf_lrr-{{ figure_baf_lrr[-1][1] }} show
{%- endif %}
the log R ratio and the B allele frequency versus the position on chromosome X
and Y for the problematic samples.
{% endif %}

{% if nb_problems > 1 %}
{{ table | safe }}

: Summarization of the gender problems encountered during Plink's analysis. HET
is the heterozygosity rate on the X chromosome. NOCALL is the percentage of no
calls on the Y chromosome. {{ "{#" }}tbl-sexcheck-results}

{% if figure_intensities %}
![
    Gender check using Plink. Mean $x$ and $y$ intensities are shown for each
    sample. Males are shown in blue, and females in red. Triangles show
    problematic samples (green for males, mauve for females). Unknown gender
    are shown in gray.
]({{ figure_intensities }}){{ "{#" }}fig-sexcheck-intensities}
{% endif %}

{% if figure_baf_lrr|length > 0 %}
{% for figure_path, sample_id in figure_baf_lrr %}
![
    Plots showing the log R ratio and the B allele frequency for chromosome X
    and Y (on the left and right, respectively) for sample {{ sample_id }}.
]({{ figure_path }}){{ "{#" }}fig-baf_lrr-{{ sample_id }}}
{% endfor %}
{% endif %}
{% endif %}
"""

        # The summary template
        self.summary_template = Environment(loader=BaseLoader)\
            .from_string(summary)

    def generate_summary(self) -> str:
        """Generate the summary from the arguments."""
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

        return self.summary_template.render(
            male_f=self.args.male_f,
            female_f=self.args.female_f,
            nb_problems=nb_problems,
            table=sex_problems.to_markdown(index=False),
            figure_intensities=figure_intensities,
            figure_baf_lrr=list(zip(baf_lrr_figures, baf_lrr_samples)),
        )
