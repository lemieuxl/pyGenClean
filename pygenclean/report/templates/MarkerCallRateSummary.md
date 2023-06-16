{% extends "section_template.md" %}

{% block section_content %}

Using a `geno` threshold of {{ geno }} (_i.e._ keeping only markers with a
missing rate $\leq{{ geno }}$), {{ "{:,d}".format(nb_markers) }}
marker{{ "s" if nb_markers > 1 }} {{ "was" if nb_markers == 1 else "were" }}
excluded from the dataset.

{% endblock %}
