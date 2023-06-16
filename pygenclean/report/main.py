"""Main report."""


import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Union

from jinja2 import BaseLoader, Environment

from ..pipeline.tree import Tree
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

```{dot}
//| fig-cap: Pipeline summary.
//| label: fig-pipeline-summary
//| fig-width: 100%
digraph pipeline {
    graph [rankdir=TB]
    node [fontsize=10 shape="box"];
    edge [fontsize=10 arrowsize=0.5];

    {% for node in dot_nodes %}
    {{ node }}
    {%- endfor %}

    {% for edge in dot_edges %}
    {{ edge }}
    {%- endfor %}
}
```

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
        dot_nodes=_generate_dot_nodes(
            initial_dataset={
                "nb_markers": nb_markers,
                "nb_samples": nb_samples,
            },
            final_datasets=kwargs["final_datasets"],
        ),
        dot_edges=_generate_dot_edges(
            tree=kwargs["qc_tree"],
            final_datasets=kwargs["final_datasets"].keys(),
        ),
    )


def _generate_dot_nodes(
    initial_dataset: Dict[str, int],
    final_datasets: Dict[str, Dict[str, str]],
) -> List[str]:
    """Generate the DOT nodes."""
    nodes = []

    # Initial dataset
    nodes.append(
        f'0 [label="{initial_dataset["nb_samples"]:,d} samples\\n'
        f'{initial_dataset["nb_markers"]:,d} markers" fillcolor="black" '
        f'style="filled" fontcolor="white" shape="note" color="white"];'
    )

    # Final datasets
    for i, (_, dataset_info) in enumerate(sorted(final_datasets.items(),
                                                 key=lambda x: int(x[0]))):
        # The number of samples
        with open(dataset_info["bfile"] + ".fam") as f:
            nb_samples = len(f.read().splitlines())

        # The number of markers
        with open(dataset_info["bfile"] + ".bim") as f:
            nb_markers = len(f.read().splitlines())

        # The node
        node = (
            f'FINAL{i + 1} [label="{nb_samples:,d} samples\\n{nb_markers:,d} '
            f'markers" fillcolor="black" style="filled" fontcolor="white" '
            f'shape="note" color="white"];'
        )
        nodes.append(node)

    return nodes


def _generate_dot_edges(tree: Tree, final_datasets: List) -> List[str]:
    """Generate the DOT edges."""
    # The final edges
    edges = []

    # The nodes and edges we've already precessed
    done_edges = set()
    done_nodes = set()

    # Processing all final datasets
    for i, final_dataset in enumerate(sorted(final_datasets, key=int)):
        # Adding the final dataset key
        edges.append(f'{final_dataset} -> FINAL{i + 1} [color="red"]')

        # Adding the nodes up to the root
        for node in tree.get_from_node_to_root(final_dataset):
            # Direct edges (from parent)
            if node.parent:
                edge = (node.parent, node.name)
                if edge not in done_edges:
                    done_edges.add(edge)
                    edges.append(f'{edge[0]} -> {edge[1]} [color="red"];')

            # Data from edges
            for sub_node in node.data_from:
                edge = (sub_node, node.name)
                if edge not in done_edges:
                    done_edges.add(edge)
                    edges.append(f'{edge[0]} -> {edge[1]} [style="dotted"];')

            # We're done with this node
            done_nodes.add(node.name)

    # Processing missing nodes
    for missing_step in tree.get_nodes() - done_nodes:
        for node in tree.get_from_node_to_root(missing_step):
            if node.parent:
                edge = (node.parent, node.name)
                if edge not in done_edges:
                    done_edges.add(edge)
                    edges.append(f"{edge[0]} -> {edge[1]};")

    return sorted(edges, key=_edge_sort)


_EDGE_RE = re.compile(r"(\d+) -> ((\d+)|(\S+))")


def _edge_sort(edge: str):
    """Edge sorter."""
    match = _EDGE_RE.match(edge)
    if match.group(3):
        return (int(match.group(1)), int(match.group(3)))
    return (int(match.group(1)), math.inf)
