#!/usr/bin/env python2.7

# This file is part of pyGenClean.
#
# pyGenClean is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# pyGenClean is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# pyGenClean.  If not, see <http://www.gnu.org/licenses/>.


import os
import re
from datetime import datetime

from . import utils as latex
from ..pipeline_error import ProgramError
from .. import __version__ as pygenclean_version


def create_report(outdirname, report_filename, **kwargs):
    """Creates a LaTeX report.

    :param report_filename: the name of the file.
    :param outdirname: the name of the output directory.

    :type report_filename: str
    :type outdirname: str

    """
    # Checking the required variables
    if "steps" in kwargs:
        assert "descriptions" in kwargs
        assert "long_descriptions" in kwargs
        assert "steps_filename" not in kwargs
    else:
        assert "steps_filename" in kwargs
        assert "descriptions" not in kwargs
        assert "long_descriptions" not in kwargs
    assert "summaries" in kwargs
    assert "background" in kwargs
    assert "project_name" in kwargs
    assert "summary_fn" in kwargs
    assert "report_title" in kwargs
    assert "report_author" in kwargs
    assert "initial_files" in kwargs
    assert "final_nb_markers" in kwargs
    assert "final_nb_samples" in kwargs
    assert "final_files" in kwargs
    assert "plink_version" in kwargs
    assert "graphic_paths_fn" in kwargs

    # Formatting the background section
    background_section = _format_background(kwargs["background"])

    # Writing the method steps to a separate file (for access later)
    steps_filename = None
    if "steps_filename" in kwargs:
        steps_filename = kwargs["steps_filename"]
    else:
        steps_filename = os.path.join(outdirname, "steps_summary.tex")
        with open(steps_filename, "w") as o_file:
            zipped = zip(kwargs["steps"], kwargs["descriptions"],
                         kwargs["long_descriptions"])
            for step, desc, long_desc in zipped:
                if desc.endswith("."):
                    desc = desc[:-1]
                step = step.replace("_", r"\_")
                to_print = latex.item(desc)
                to_print += " [{}].".format(latex.texttt(step))
                if long_desc is not None:
                    to_print += " " + long_desc
                print >>o_file, latex.wrap_lines(to_print) + "\n"

    # Adding the content of the results section
    result_summaries = []
    for name in kwargs["summaries"]:
        full_path = os.path.abspath(name)
        if os.path.isfile(full_path):
            rel_path = os.path.relpath(full_path, outdirname)
            result_summaries.append(re.sub(r"\\", "/", rel_path))

    # Reading the initial_files file
    initial_files = None
    with open(kwargs["initial_files"], "r") as i_file:
        initial_files = i_file.read().splitlines()

    # Reading the final_files file
    final_files = None
    with open(kwargs["final_files"], "r") as i_file:
        final_files = [i.split("\t")[0] for i in i_file.read().splitlines()]

    # Adding the bibliography content
    biblio_entry = latex.bib_entry(
        name="pyGenClean",
        authors="Lemieux Perreault LP, Provost S, Legault MA, Barhdadi A, "
                r"Dub\'e MP",
        title="pyGenClean: efficient tool for genetic data clean up before "
              "association testing",
        journal="Bioinformatics",
        year="2013",
        volume="29",
        number="13",
        pages="1704--1705",
    ) + "\n" * 2 + latex.bib_entry(
        name="plink",
        authors="Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, "
                "Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ, Sham PC",
        title="PLINK: a tool set for whole-genome association and "
              "population-based linkage analyses",
        journal="American Journal of Human Genetics",
        year="2007",
        volume="81",
        number="3",
        pages="559--575",
    ) + "\n" * 2 + latex.bib_entry(
        name="bafRegress",
        authors=r"Goo J, Matthew F, Kurt NH, Jane MR, Kimberly FD, "
                r"Gon{\c{c}}alo RA, Michael B, Hyun Min K",
        title="Detecting and estimating contamination of human DNA samples in "
              "sequencing and array-based genotype data",
        journal="The American Journal of Human Genetics",
        year="2012",
        volume="91",
        number="5",
        pages="839--848",
    )

    # Getting the template
    main_template = latex.jinja2_env.get_template("main_document.tex")

    # Getting the data
    today = datetime.today()

    # Reading the graphics path
    graphic_paths = []
    if kwargs["graphic_paths_fn"] is not None:
        with open(kwargs["graphic_paths_fn"], "r") as i_file:
            graphic_paths = [
                re.sub(r"\\", "/", path) + ("" if path.endswith("/") else "/")
                for path in i_file.read().splitlines()
            ]

    try:
        with open(report_filename, "w") as i_file:
            # Rendering the template
            print >>i_file, main_template.render(
                project_name=latex.sanitize_tex(kwargs["project_name"]),
                month=today.strftime("%B"),
                day=today.day,
                year=today.year,
                background_content=background_section,
                result_summaries=result_summaries,
                bibliography_content=biblio_entry,
                pygenclean_version=pygenclean_version,
                plink_version=kwargs["plink_version"],
                steps_filename=os.path.basename(steps_filename),
                final_results=_create_summary_table(
                    kwargs["summary_fn"],
                    latex.jinja2_env.get_template("summary_table.tex"),
                    nb_samples=kwargs["final_nb_samples"],
                    nb_markers=kwargs["final_nb_markers"],
                ),
                report_title=latex.sanitize_tex(kwargs["report_title"]),
                report_author=latex.sanitize_tex(kwargs["report_author"]),
                initial_files=initial_files,
                final_files=final_files,
                final_nb_samples=kwargs["final_nb_samples"],
                final_nb_markers=kwargs["final_nb_markers"],
                graphic_paths=graphic_paths,
            )

    except IOError:
        msg = "{}: could not create report".format(report_filename)
        raise ProgramError(msg)


def _format_background(background):
    """Formats the background section

    :param background: the background content or file.

    :type background: str or file

    :returns: the background content.
    :rtype: str

    """
    # Getting the background
    if os.path.isfile(background):
        with open(background, "r") as i_file:
            background = i_file.read().splitlines()
    else:
        background = background.splitlines()

    # Formatting
    final_background = ""
    for line in background:
        if line == "":
            final_background += r"\\" + "\n\n"
            continue

        final_background += latex.wrap_lines(latex.sanitize_tex(line))

    return final_background


def _create_summary_table(fn, template, nb_samples, nb_markers):
    """Creates the final table.

    :param fn: the name of the file containing the summary.
    :param template: the Jinja2 template.
    :param nb_samples: the final number of samples.
    :param nb_markers: the final number of markers.

    :type fn: str
    :type template: Jinja2.template
    :type nb_samples: str
    :type nb_markers: str

    """
    # The final data
    table_data = []

    # Reading the summary file
    with open(fn, "r") as i_file:
        data = None

        line = i_file.readline()
        while line != "":
            if line.startswith("#"):
                # If there is data, this isn't the first line, so we save
                if data:
                    table_data.append(data)

                # This is the 'header' of a section (hence a new section)
                data = dict(
                    header=line.rstrip("\r\n").split(" ")[1],
                    data=[],
                )

                # Changing to next line
                line = i_file.readline()
                continue

            # If the line starts with '---', then it's a horizontal line
            if line.startswith("---"):
                data["data"].append(dict(hline=True))

                # Changing to next line
                line = i_file.readline()
                continue

            # If the line starts with '  -', then it's a sub section
            if line.startswith("  -"):
                tmp = line[4:].rstrip("\r\n").split("\t")
                if data["header"].endswith("/subset"):
                    if tmp[0].startswith("_file_path:"):
                        tmp[0] = r"\path{" + tmp[0][11:] + "}"
                elif data["header"].endswith("/flag_hw"):
                    tmp[0] = latex.format_numbers(tmp[0], prefix="p < ")
                else:
                    tmp = map(latex.sanitize_tex, tmp)
                    if tmp[0].startswith("x"):
                        tmp[0] = latex.inline_math(r"\times " + tmp[0][1:])

                data["data"].append(dict(
                    hline=False,
                    multicol=False,
                    row_data=tmp,
                ))

                # Changing to next line
                line = i_file.readline()
                continue

            # This is a regular line
            data["data"].append(dict(
                hline=False,
                multicol=True,
                row_data=map(
                    latex.sanitize_tex,
                    line.rstrip("\r\n").split("\t"),
                ),
            ))

            # Skipping to next line
            line = i_file.readline()

    # We add the last entry
    table_data.append(data)

    # Rendering
    return template.render(table_data=table_data, final_nb_markers=nb_markers,
                           final_nb_samples=nb_samples)
