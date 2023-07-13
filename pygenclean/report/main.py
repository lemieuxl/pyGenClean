"""Main report."""


import logging
import math
import re
from pathlib import Path
from typing import IO, Dict, List, Optional, Tuple, Union

from jinja2 import BaseLoader, Environment
from tabulate import tabulate

from ..pipeline.tree import QCNode, Tree
from ..qc_modules import qc_module_impact
from ..qc_modules.marker_call_rate.marker_call_rate import _DEFAULT_GENO
from ..qc_modules.sample_call_rate.sample_call_rate import _DEFAULT_MIND
from ..utils import count_lines
from ..utils.plink import compare_bim, compare_fam
from ..utils.task import execute_external_command
from ..version import pygenclean_version


logger = logging.getLogger(__name__)


MAIN_TEMPLATE = Environment(loader=BaseLoader).from_string("""\
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

@fig-pipeline-summary-dot shows the summary of the pipeline, and
@tbl-final-summary-general shows the general steps' summary.

{% if is_dot_code %}
```{dot}
//| fig-cap: Pipeline summary.
//| label: fig-pipeline-summary-dot
//| fig-width: 100%
{{ pipeline_summary|safe }}
```
{% else %}
![
    Pipeline summary.
]({{ pipeline_summary }}){{ "{#" }}fig-pipeline-summary-dot}
{% endif %}

{{ conclusion_summaries["other_steps"]["summary_table"]|safe }}

: Summary information of the data cleanup procedures for the general steps.
{{ "{#" }}tbl-final-summary-general}

{% for step, step_info in final_datasets %}
### {{ step_info["desc"] if step_info["desc"] else "Step " + step }}

After the genetic data clean up procedure, a total of
{{ "{:,d}".format(conclusion_summaries[step]["nb_samples"]) }} samples and
{{ "{:,d}".format(conclusion_summaries[step]["nb_markers"]) }} markers
remained. The following files are available for downstream analysis, and
@tbl-final-summary-{{ step }} shows a summary of exclusions for this
dataset.

- `{{ conclusion_summaries[step]["prefix"] }}.bed`
- `{{ conclusion_summaries[step]["prefix"] }}.bim`
- `{{ conclusion_summaries[step]["prefix"] }}.fam`

{{ conclusion_summaries[step]["summary_table"]|safe }}

: Summary information of the data cleanup procedure for
{{ step_info["desc"] if step_info["desc"] else "Step " + step }}.
{{ "{#" }}tbl-final-summary-{{ step }}}
{% endfor %}
""")


DOT_TEMPLATE = Environment(loader=BaseLoader).from_string("""\
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
""")


def generate_report(**kwargs: Dict[str, Optional[Union[str, int]]]) -> str:
    """Generate the report."""
    # The QC directory
    qc_dir = Path(kwargs["qc_dir"])

    # We want to compute the n's for each of the dataset
    step_stats = _compute_step_stats(
        qc_conf=kwargs["qc_conf"],
        qc_dir=qc_dir,
        tree=kwargs["qc_tree"],
        datasets=kwargs["final_datasets"],
    )

    # The number of markers and samples of the initial files
    nb_markers = count_lines(kwargs["bfile"] + ".bim")
    nb_samples = count_lines(kwargs["bfile"] + ".fam")

    # Writing the files containing the methods for each step and the results
    logger.info("Generating the methods and results sections")
    with open(qc_dir / "methods.qmd", "w") as methods_f, \
         open(qc_dir / "results.qmd", "w") as results_f:
        # Skipping first step, because it's step 0 (initial dataaset)
        for i, step_name in enumerate(kwargs["qc_tree"].get_node_order()[1:]):
            # The QC node
            qc_node = kwargs["qc_tree"].get_node(step_name)

            # The methods for this step
            qc_methods = qc_node.summary.generate_methods()
            print(f"{i + 1}. {qc_methods} [`{qc_node.module_name}`]",
                  file=methods_f)

            # The results for this step
            filename = qc_node.summary.write_results()
            print("{{< include " + str(filename) + " >}}\n", file=results_f)

    # Generating the pipeline summary figure
    pipeline_summary = _generate_summary_figure(
        qc_dir=qc_dir,
        dot_nodes=_generate_dot_nodes(
            initial_dataset={
                "nb_markers": nb_markers,
                "nb_samples": nb_samples,
            },
            final_datasets=kwargs["final_datasets"],
            all_steps=kwargs["qc_tree"].get_nodes(),
            step_stats=step_stats,
            step_conf=kwargs["qc_conf"]["steps"],
        ),
        dot_edges=_generate_dot_edges(
            tree=kwargs["qc_tree"],
            final_datasets=kwargs["final_datasets"].keys(),
        ),
        use_dot=kwargs["use_dot"],
    )

    return MAIN_TEMPLATE.render(
        docx_template=kwargs["report_template"],
        title=kwargs["report_title"],
        subtitle=kwargs["report_number"],
        authors=", ".join(kwargs["report_authors"]),
        version=pygenclean_version,
        bfile=kwargs["bfile"],
        nb_markers=nb_markers,
        nb_samples=nb_samples,
        step_methods=str(qc_dir / "methods.qmd"),
        step_results=str(qc_dir / "results.qmd"),
        pipeline_summary=pipeline_summary,
        is_dot_code=not kwargs["use_dot"],
        final_datasets=sorted(kwargs["final_datasets"].items(),
                              key=lambda x: int(x[0])),
        conclusion_summaries=_generate_conlusion_summaries(
            initial_numbers={
                "nb_markers": nb_markers,
                "nb_samples": nb_samples,
            },
            tree=kwargs["qc_tree"],
            datasets=kwargs["final_datasets"],
            stats=step_stats,
            conf=kwargs["qc_conf"]["steps"],
        ),
    )


def _generate_conlusion_summaries(
    initial_numbers: Dict[str, int],
    tree: Tree,
    datasets: Dict[str, Dict[str, str]],
    stats: Dict[str, Dict[str, int]],
    conf: dict,
) -> Dict[str, Dict[str, str]]:
    """Generate the final dataset summary (Markdown tables)"""
    conclusion_summaries = {}

    summarized_steps = set()

    for dataset in datasets.keys():
        dataset_summary = []

        for qc_node in reversed(tree.get_from_node_to_root(dataset)):
            if qc_node.name == "0":
                assert len(dataset_summary) == 0

                # Initial numbers
                dataset_summary.append((
                    None,
                    "Initial numbers",
                    None,
                    initial_numbers['nb_markers'],
                    initial_numbers['nb_samples'],
                ))
                continue

            step_info = stats[qc_node.name]

            qc_module = conf[qc_node.name]["module"].replace("-", "_")

            description = qc_node.summary.result_section_name

            # Is this a subset?
            if qc_module == "subset":
                description += " - " + conf[qc_node.name]["reason"]

            # Is this a sample call rate?
            elif qc_module == "sample_call_rate":
                description += (
                    " - mind="
                    + str(conf[qc_node.name].get("mind", _DEFAULT_MIND))
                )

            # Is this a marker call rate?
            elif qc_module == "marker_call_rate":
                description += (
                    " - geno="
                    + str(conf[qc_node.name].get("geno", _DEFAULT_GENO))
                )

            summary_table_info = qc_node.summary.summary_table_info
            if summary_table_info:
                assert len(summary_table_info) == 1
                summary_table_info = summary_table_info[0][-1]

            dataset_summary.append((
                int(qc_node.name),
                description,
                summary_table_info,
                -step_info["nb_markers"] if step_info["nb_markers"] else None,
                -step_info["nb_samples"] if step_info["nb_samples"] else None,
            ))

            # This step was summarized
            summarized_steps.add(qc_node.name)

            # Is this the final step?
            if dataset == qc_node.name:
                dataset_summary.append((
                    None,
                    "Final numbers",
                    None,
                    stats[dataset]["bim"],
                    stats[dataset]["fam"],
                ))

        conclusion_summaries[dataset] = {
            "summary_table": tabulate(
                dataset_summary,
                headers=("Step", "Description", "Info", "Markers", "Samples"),
                intfmt=",",
                tablefmt="github",
            ),
            "nb_markers": dataset_summary[-1][-2],
            "nb_samples": dataset_summary[-1][-1],
            "prefix": tree.get_node(dataset).bfile,
        }

    # Adding the additionl summaries which are not present in the data
    dataset_summary = []
    for step_name in tree.get_node_order():
        if step_name == "0" or step_name in summarized_steps:
            continue

        qc_node = tree.get_node(step_name)

        # Generating the description and the information
        description = qc_node.summary.result_section_name
        information = None

        if qc_node.summary.summary_table_info:
            information = ""
            for i, (desc, info) in enumerate(
                qc_node.summary.summary_table_info
            ):
                if i == 0:
                    description += "\n"
                description += "\n  - " + desc
                information += f"\n- {info:,d}"

        dataset_summary.append((
            int(qc_node.name),
            int(qc_node.parent),
            description,
            information,
        ))

    conclusion_summaries["other_steps"] = {
        "summary_table": tabulate(
            dataset_summary,
            headers=("Step", "Parent Step", "Description", "Information"),
            intfmt=",",
            tablefmt="grid",
        ),
    }

    return conclusion_summaries


def _generate_summary_figure(qc_dir: Path, dot_nodes: List[str],
                             dot_edges: List[str], use_dot: bool) -> str:
    """Generate the summary figure."""
    dot_code = DOT_TEMPLATE.render(dot_nodes=dot_nodes, dot_edges=dot_edges)

    # If we don't have to use dot, we just return the dot code (quarto will do
    # the rest)
    if not use_dot:
        return dot_code

    # We need to generate the dot file and compile it using dot
    dot_file = qc_dir / "pipeline_summary.dot"
    with open(dot_file, "w") as f:
        print(dot_code, file=f)

    # We execute dot
    png_file = qc_dir / "pipeline_summary.png"
    with open(png_file, "wb") as f:
        f.write(execute_external_command(["dot", "-Tpng", str(dot_file)],
                                         decode=False))

    return str(png_file)


def _compute_step_stats(
        qc_dir: Path,
        qc_conf: dict,
        tree: Tree,
        datasets: Dict[str, Dict[str, str]],
) -> Dict[str, Dict[str, int]]:
    """Compute datasets statistics"""
    step_stats = {}

    for dataset, dataset_info in sorted(datasets.items(),
                                        key=lambda x: int(x[0])):
        logger.info("Generating exclusion lists for dataset %s (%s)",
                    dataset, dataset_info["desc"])
        # The two files (excluded markers and samples)
        excluded_markers = qc_dir / f"excluded_markers_step_{dataset}.txt"
        excluded_samples = qc_dir / f"excluded_samples_step_{dataset}.txt"

        total_nb_samples = 0
        total_nb_markers = 0
        with open(excluded_markers, "w") as f_excluded_markers,\
             open(excluded_samples, "w") as f_excluded_samples:
            previous_qc_node = None
            for qc_node in reversed(tree.get_from_node_to_root(dataset)):
                if not previous_qc_node:
                    previous_qc_node = qc_node
                    continue

                # Comparing previous (newest) with current step
                # Comparing FAM
                nb_samples, nb_samples_bfile = _write_sample_exclusions(
                    qc_node=qc_node,
                    previous_qc_node=previous_qc_node,
                    step_conf=qc_conf["steps"][qc_node.name],
                    f=f_excluded_samples,
                )

                # Comparing BIM
                nb_markers, nb_markers_bfile = _write_marker_exclusions(
                    qc_node=qc_node,
                    previous_qc_node=previous_qc_node,
                    step_conf=qc_conf["steps"][qc_node.name],
                    f=f_excluded_markers,
                )

                total_nb_samples += nb_samples
                total_nb_markers += nb_markers

                if qc_node.name in step_stats:
                    assert step_stats[qc_node.name]["nb_samples"] == nb_samples
                    assert step_stats[qc_node.name]["nb_markers"] == nb_markers
                    assert step_stats[qc_node.name]["fam"] == nb_samples_bfile
                    assert step_stats[qc_node.name]["bim"] == nb_markers_bfile

                else:
                    step_stats[qc_node.name] = {
                        "nb_samples": nb_samples,
                        "nb_markers": nb_markers,
                        "fam": nb_samples_bfile,
                        "bim": nb_markers_bfile
                    }

                # Next step
                previous_qc_node = qc_node

        logger.info("  - %d markers were excluded", total_nb_markers)
        logger.info("  - %d samples were excluded", total_nb_samples)

    return step_stats


def _write_marker_exclusions(qc_node: QCNode, previous_qc_node: QCNode,
                             step_conf: str, f: IO) -> Tuple[int, int]:
    """Write the marker exclusion to file."""
    before_only, both, after_only = compare_bim(
        previous_qc_node.bfile + ".bim", qc_node.bfile + ".bim",
    )

    # Making sure no markers were added...
    assert len(after_only) == 0

    if before_only:
        # Some markers were excluded
        exclusion_info = _get_exclusion_info(
            step=qc_node.name,
            qc_module=step_conf["module"],
            reason=step_conf.get("reason"),
        )
        for marker in sorted(before_only):
            print(marker, exclusion_info, sep="\t", file=f)

    return len(before_only), len(both)


def _write_sample_exclusions(qc_node: QCNode, previous_qc_node: QCNode,
                             step_conf: str, f: IO) -> Tuple[int, int]:
    """Write the sample exclusion to file."""
    before_only, both, after_only = compare_fam(
        previous_qc_node.bfile + ".fam", qc_node.bfile + ".fam",
    )

    # Making sure no samples were added...
    assert len(after_only) == 0

    if before_only:
        # Some samples were excluded
        exclusion_info = _get_exclusion_info(
            step=qc_node.name,
            qc_module=step_conf["module"],
            reason=step_conf.get("reason"),
        )
        for fid, iid in sorted(before_only):
            print(fid, iid, exclusion_info, sep="\t", file=f)

    return len(before_only), len(both)


def _get_exclusion_info(step: str, qc_module: str,
                        reason: Optional[str]) -> str:
    """Get the exclusion information."""
    exclusion_info = f"{step} {qc_module.replace('-', '_')}"
    if reason:
        exclusion_info += f" ({reason})"
    return exclusion_info


def _generate_dot_nodes(
    initial_dataset: Dict[str, int],
    final_datasets: Dict[str, Dict[str, str]],
    all_steps: List[str],
    step_stats: Dict[str, Dict[str, int]],
    step_conf: dict,
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
    for i, (step, dataset_info) in enumerate(sorted(final_datasets.items(),
                                                    key=lambda x: int(x[0]))):
        # The description of the dataset
        description = dataset_info["desc"]
        if description:
            description += r"\n"
        else:
            description = ""

        # The number of samples and markers
        nb_samples = step_stats[step]["fam"]
        nb_markers = step_stats[step]["bim"]

        nodes.append(
            f'FINAL{i + 1} [label="{description}{nb_samples:,d} samples\\n'
            f'{nb_markers:,d} markers" fillcolor="black" style="filled" '
            f'fontcolor="white" shape="note" color="white"];'
        )

    for step in all_steps:
        if step == "0":
            continue

        # The name of the QC module
        qc_module = step_conf[step]["module"].replace("-", "_")

        # The node label
        node_label = _create_node_label(
            step=step,
            step_conf=step_conf[step],
            step_stats=step_stats.get(step),
        )

        # The node style
        node_style = _create_node_style(step_stats.get(step))

        # The node color (light gray when impact on markers)
        color = _create_node_color(qc_module, node_label)

        nodes.append(
            f'{step} [label="{node_label}"{color}{node_style}];'
        )

    return nodes


def _create_node_style(stats: Optional[Dict[str, int]]) -> str:
    """Create the node style."""
    if stats:
        if sum(stats.values()) > 0:
            return ""
    return ' style="rounded"'


_SUBSET_MARKER_RE = re.compile(r"\\n-[0-9,]+ markers$")


def _create_node_color(qc_module: str, label: str) -> str:
    """Create the node color."""
    if (
        (qc_module_impact[qc_module] == "markers")
        or (qc_module == "subset" and _SUBSET_MARKER_RE.search(label))
    ):
        return ' fillcolor="lightgray" style="filled"'

    return ""


def _create_node_label(step: str, step_conf: dict,
                       step_stats: Optional[Dict[str, int]]) -> str:
    """create the node label."""
    # The QC module
    qc_module = step_conf["module"].replace("-", "_")

    # Reason if any
    reason = step_conf.get("reason", "")
    if reason:
        reason = r"\n" + reason

    # Parameters if marker_call_rate sample_call_rate
    parameter = ""
    if qc_module == "marker_call_rate":
        parameter = r"\ngeno=" + str(step_conf.get("geno", _DEFAULT_GENO))
    elif qc_module == "sample_call_rate":
        parameter = r"\nmind=" + str(step_conf.get("mind", _DEFAULT_MIND))

    # Step diff (markers or samples)
    step_diff = []
    for item_type in ("markers", "samples"):
        if step_stats:
            nb_item = step_stats["nb_" + item_type]
            if nb_item:
                step_diff.append(f"-{nb_item:,d} {item_type}")
    step_diff = ", ".join(step_diff)
    if step_diff:
        step_diff = r"\n" + step_diff

    return f"{step}. {qc_module}{reason}{parameter}{step_diff}"


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
