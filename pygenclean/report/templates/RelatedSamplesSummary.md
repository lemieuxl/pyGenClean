{% extends "section_template.md" %}

{% block section_content %}

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
@fig-{{ label_prefix }}-z1 shows $Z_1$ versus $IBS2_{ratio}^\ast$ for all
related samples found by _Plink_.
{% endif -%}
{% if figure_z2 -%}
@fig-{{ label_prefix }}-z2 shows $Z_2$ versus $IBS2_{ratio}^\ast$ for all
related samples found by _Plink_.
{% endif -%}
{% if table -%}
@tbl-{{ label_prefix }}-merged-related-samples lists the related sample pairs
with estimated relationship.
{% endif %}

{% if figure_z1 %}
![
    $Z_1$ versus $IBS2_{ratio}^\ast$ for all related samples found by _Plink_
    in the IBS analysis.
]({{ figure_z1 }}){{ "{#" }}fig-{{ label_prefix }}-z1}
{% endif %}

{% if figure_z2 %}
![
    $Z_2$ versus $IBS2_{ratio}^\ast$ for all related samples found by _Plink_
    in the IBS analysis.
]({{ figure_z2 }}){{ "{#" }}fig-{{ label_prefix }}-z2}
{% endif %}

{% if table %}
{{ table }}

: List of all related samples with estimated relationship. Sample pairs are
grouped according to their estimated family (the index
column). {{ "{#" }}tbl-{{ label_prefix }}-merged-related-samples}
{% endif %}

{%- endif %}

{% endblock %}
