{% extends "section_template.md" %}

{% block section_content %}

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

{% endblock %}
