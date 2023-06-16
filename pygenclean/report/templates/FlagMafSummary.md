{% extends "section_template.md" %}

{% block section_content %}

After computing minor allele frequencies (MAF) of all markers using _Plink_, a
total of {{ "{:,d}".format(nb_flagged) }} marker{{ "s" if nb_flagged > 1 }} had
a MAF of zero and were flagged (see file `flag_maf_0.list` for more
information).

{% endblock %}
