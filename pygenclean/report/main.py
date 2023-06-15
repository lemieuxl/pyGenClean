"""Main report."""


from pathlib import Path
from typing import Dict, Optional, Union

from jinja2 import BaseLoader, Environment

from ..version import pygenclean_version


TEMPLATE = Environment(loader=BaseLoader).from_string("""\
---
title: "{{ title }}"
subtitle: "{{ subtitle }}"
author: "{{ authors }}"
date: today
date-format: long
format:
    docx:
        toc: true
        number-sections: true
        {% if docx_template -%}
        reference-doc: {{ docx_template }}
        {%- endif %}
---

## Background

## Methods

The automated script _pyGenClean_ version {{ version }} was used to launch the
analysis. Plink was used for the data cleanup procedure. The following files
were used as input. They contained {{ "{:,d}".format(nb_samples) }} samples
genotyped for {{ "{:,d}".format(nb_markers) }} genetic markers.

- `{{ bfile }}.bed`
- `{{ bfile }}.bim`
- `{{ bfile }}.fam`

Briefly, the clean up procedure was as follow:

{{ "{{< include " + step_methods + " >}}" }}

## Results

{{ "{{< include " + step_results + " >}}" }}

## Conclusions

After the genetic data clean up procedure, a total of
{{ "{:,d}".format(final_nb_samples) }} samples and
{{ "{:,d}".format(final_nb_markers) }} markers remained. The following files
are available for downstream analysis.

- `{{ final_prefix }}.bed`
- `{{ final_prefix }}.bim`
- `{{ final_prefix }}.fam`
""")


def generate_report(**kwargs: Dict[str, Optional[Union[str, int]]]) -> str:
    """Generate the report."""
    # The QC directoru
    qc_dir = Path(kwargs["qc_dir"])

    # The number of markers
    with open(kwargs["bfile"] + ".bim") as f:
        nb_markers = len(f.read().splitlines())

    # The number of samples
    with open(kwargs["bfile"] + ".fam") as f:
        nb_samples = len(f.read().splitlines())

    # Writing the file containing the methods for each step
    methods_file = qc_dir / "methods.qmd"
    with open(methods_file, "w") as f:
        zipped = zip(kwargs["qc_modules"], kwargs["step_methods"])
        for i, (qc_module, qc_methods) in enumerate(zipped):
            print(f"{i + 1}. {qc_methods} [`{qc_module}`]", file=f)

    # Writing the file containing the results for each step
    results_file = qc_dir / "results.qmd"
    with open(results_file, "w") as f:
        for file_name in kwargs["step_results"]:
            print("{{< include " + str(file_name) + " >}}\n", file=f)

    # Getting the number of samples and markers from the final file
    with open(kwargs["final_prefix"] + ".bim") as f:
        final_nb_markers = len(f.read().splitlines())
    with open(kwargs["final_prefix"] + ".fam") as f:
        final_nb_samples = len(f.read().splitlines())

    return TEMPLATE.render(
        docx_template=kwargs["report_template"],
        title=kwargs["report_title"],
        subtitle=kwargs["report_number"],
        authors=", ".join(kwargs["report_authors"]),
        version=pygenclean_version,
        bfile=kwargs["bfile"],
        nb_markers=nb_markers,
        nb_samples=nb_samples,
        final_nb_markers=final_nb_markers,
        final_nb_samples=final_nb_samples,
        final_prefix=kwargs["final_prefix"],
        step_methods=str(methods_file),
        step_results=str(results_file),
    )
