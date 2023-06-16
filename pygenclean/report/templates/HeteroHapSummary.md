{% extends "section_template.md" %}

{% block section_content %}

After _Plink_'s heterozygous haploid analysis, a total of
{{ "{:,d}".format(nb_hetero_hap) }} genotype{{ "s" if nb_hetero_hap > 1}}
{{ "was" if nb_hetero_hap == 1 else "were" }} set to missing.

{% endblock %}
