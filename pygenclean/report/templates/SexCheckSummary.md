{% extends "section_template.md" %}

{% block section_content %}

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

: Summarization of the sex mismatch problems encountered during Plink's
analysis. `HET` is the heterozygosity rate on the X chromosome. `NOCALL` is the
percentage of no calls on the Y
chromosome. {{ "{#" }}tbl-{{ label_prefix }}-results}

{% if figure_intensities %}
![
    Sex check using Plink. Mean $x$ and $y$ intensities are shown for each
    sample. Males are shown in blue, and females in red. Triangles show
    problematic samples (green for males, mauve for females). Unknown sex
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

{% endblock %}
