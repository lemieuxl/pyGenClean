{% extends "section_template.md" %}

{% block section_content %}

After scrutiny, {{ "{:,d}".format(all_failed) }}
marker{{ "s" if all_failed > 1 }} {{ "was" if all_failed == 1 else "were" }}
excluded from the dataset because of a call rate of 0.
{%- if not keep_het %}
Also, {{ "{:,d}".format(all_hetero) }} marker{{ "s" if all_hetero > 1 }}
{{ "was" if all_hetero == 1 else "were" }} excluded from the dataset because all
samples were heterozygous (excluding the mitochondrial chromosome).
{%- else %}
Also, {{ "{:,d}".format(all_hetero) }} marker{{ "s" if all_hetero > 1 }}
{{ "was" if all_hetero == 1 else "were" }} found to only have heterozygous
genotypes, but they were kept in the dataset (because of the `--keep-het`
flag).
{% endif %}

{% endblock %}
