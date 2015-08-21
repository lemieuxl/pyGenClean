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
from datetime import datetime

import pyGenClean.LaTeX as latex
from pyGenClean import __version__ as pygenclean_version


def create_report(outdirname, report_filename, **kwargs):
    """Creates a LaTeX report.

    :param report_filename: the name of the file.
    :param outdirname: the name of the output directory.

    :type report_filename: string
    :type outdirname: string

    """
    # Checking the required variables
    assert "steps" in kwargs
    assert "logo_path" in kwargs
    assert "summaries" in kwargs
    assert "background" in kwargs
    assert "descriptions" in kwargs
    assert "project_name" in kwargs

    # Getting the required information
    project_name = kwargs["project_name"]
    logo_path = kwargs["logo_path"]
    prog_version = ".".join(pygenclean_version)
    steps = kwargs["steps"]
    summaries = kwargs["summaries"]
    descriptions = kwargs["descriptions"]
    background_section = latex.wrap_lines(kwargs["background"])

    # Writing the method steps to a separate file (for access later)
    step_filename = os.path.join(outdirname, "steps_summary.tex")
    with open(step_filename, "w") as o_file:
        for step, desc in zip(steps, descriptions):
            step = step.replace("_", r"\_")
            to_print = latex.item(desc)
            to_print += " ({})".format(latex.texttt(step))
            print >>o_file, latex.wrap_lines(to_print) + "\n"

    # Adding the content of the results section
    result_summaries = []
    for name in summaries:
        full_path = os.path.abspath(name)
        if os.path.isfile(full_path):
            rel_path = os.path.relpath(full_path, outdirname)
            result_summaries.append(rel_path)

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
    )

    # Getting the template
    main_template = latex.jinja2_env.get_template("main_document.tex")

    # Getting the data
    today = datetime.today()

    try:
        with open(report_filename, "w") as i_file:
            # Rendering the template
            # The report information
            #   - project_name: the name of the project
            #   - logo_path: the path to the logo
            #   - background_content: the content of the 'background'
            #   - methods_content: the content of the 'methods'
            #   - results_content: the content of the 'results'
            #   - conclusions_content: the content of the 'conclusions'
            #   - bibliography_content: the content of the 'bibliography'
            print >>i_file, main_template.render(
                project_name=project_name,
                month=today.strftime("%B"),
                day=today.day,
                year=today.year,
                logo_dir=os.path.abspath(os.path.dirname(logo_path)) + "/",
                logo_path=os.path.basename(logo_path),
                background_content=background_section,
                result_summaries=result_summaries,
                conclusions_content="",
                bibliography_content=biblio_entry,
                pygenclean_version=pygenclean_version,
                step_filename=os.path.basename(step_filename),
            )

    except IOError:
        msg = "{}: could not create report".format(report_filename)
        raise ProgramError(msg)
