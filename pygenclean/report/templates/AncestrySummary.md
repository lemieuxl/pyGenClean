{% extends "section_template.md" %}

{% set scree_text = "@fig-" + label_prefix + "-scree shows the scree plot for the principal components of the MDS analysis." %}

{% block section_content %}

{% if skip_ref_pops -%}

Principal components analysis was performed using
{{ "{:,d}".format(nb_markers) }} marker{{ "s" if nb_markers > 1 }} on the study
dataset only.
{% if scree_figure -%}
{{ scree_text }}
{% endif %}

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
{{ scree_text }}
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

{%- endif %}

{% if scree_figure %}
![
    Scree plot for the principal components of the MDS analysis.
]({{ scree_figure }}){{ "{#" }}fig-{{ label_prefix }}-scree width='10cm'}
{% endif %}

{% endblock %}
