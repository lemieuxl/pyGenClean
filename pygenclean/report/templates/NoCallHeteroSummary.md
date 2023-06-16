{% extends "section_template.md" %}

{% block section_content %}

After scrutiny, {{ "{:,d}".format(all_failed) }}
marker{{ "s" if all_failed > 1 }} {{ "was" if all_failed == 1 else "were" }}
excluded from the dataset because of a call rate of 0. Also,
{{ "{:,d}".format(all_hetero) }} marker{{ "s" if all_hetero > 1 }}
{{ "was" if all_hetero == 1 else "was" }} excluded from the dataset because all
samples were heterozygous (excluding the mitochondrial chromosome).

{% endblock %}
