{% extends "section_template.md" %}

{% block section_content %}

Using a `mind` threshold of {{ mind }} (_i.e._ keeping only samples with a
missing rate $\leq{{ mind }}$), {{ "{:,d}".format(nb_samples) }}
sample{{ "s" if nb_samples > 1 }} {{ "was" if nb_samples == 1 else "were" }}
excluded from the dataset.
{% if missing_rate_table -%}
@tbl-{{ label_prefix }}-missing-rate shows the missing rate for excluded
samples.
{%- endif %}

{% if missing_rate_table %}
{{ missing_rate_table | safe }}

: Excluded samples' missing rate. {{ "{#" }}tbl-{{ label_prefix }}-missing-rate}
{% endif %}

{% endblock %}
