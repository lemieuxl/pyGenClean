{% extends "section_template.md" %}

{% block section_content %}

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

{% endblock %}
