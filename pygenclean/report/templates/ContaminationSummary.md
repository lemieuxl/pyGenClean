{% extends "section_template.md" %}

{% block section_content %}

A total of {{ "{:,d}".format(nb_samples) }} sample{{ "s" if nb_samples > 1 }}
{{ "were" if nb_samples > 1 else "was" }} analyzed for contamination using
_bafRegress_. The analysis was performed using
{{ "{:,d}".format(nb_autosomal) }} autosomal
marker{{ "s" if nb_autosomal > 1 }}. Using a threshold of {{ threshold }},
{{ "{:,d}".format(nb_contaminated) }} sample{{ "s" if nb_contaminated > 1 }}
{{ "were" if nb_contaminated > 1 else "was" }} estimated to be contaminated.
{%- if nb_contaminated > 0 %}
@tbl-{{ label_prefix }}-results lists all the samples that were estimated to be
contaminated (_i.e._ with an estimate $>{{ threshold }}$)

{{ table }}

: List of all possible contaminated samples (_i.e._ with an estimate computed
by _bafRegress_ $>{{ threshold }}$). {{ "{#" }}tbl-{{ label_prefix }}-results}
{% endif -%}

{% endblock %}
