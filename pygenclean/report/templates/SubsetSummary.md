{% extends "section_template.md" %}

{% block section_content %}

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

{% endblock %}
