{% extends "section_template.md" %}

{% block section_content %}

A total of {{ "{:,d}".format(nb_dup_samples) }} duplicated
sample{{ "s" if nb_dup_samples > 1 }}
{{ "was" if nb_dup_samples == 1 else "were" }} found.

{% endblock %}
