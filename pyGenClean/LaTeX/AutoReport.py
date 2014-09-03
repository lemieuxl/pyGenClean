#!/usr/bin/env python2.7
## This file is part of pyGenClean.
## 
## pyGenClean is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.
## 
## pyGenClean is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License along with
## pyGenClean.  If not, see <http://www.gnu.org/licenses/>.

import os
import textwrap

import pyGenClean
import pyGenClean.LaTeX as template


def create_report(filename, **kwargs):
    """Creates a LaTeX report.
    
    :param filename: the name of the file.
    
    :type filename: string
    
    """
    # Getting the required information
    title = kwargs["title"]
    logo_path = kwargs["logo_path"]
    prog_version = ".".join(pyGenClean.get_version())
    steps = kwargs["steps"]
    summaries = kwargs["summaries"]
    descriptions = kwargs["descriptions"]

    try:
        with open(filename, "w") as i_file:
            # Writing the preamble
            print >>i_file, template.preamble.substitute(
                        project_name=title,
                        logo_path=logo_path)

            # Writing the background sections
            print >>i_file, template.section("Background")
            text = ("Lorem ipsum dolor sit amet, consectetur adipiscing elit. "
                    "Suspendisse lectus ligula, volutpat eget convallis a, "
                    "porttitor vitae est. Pellentesque ornare ipsum vitae odio "
                    "sodales, eu elementum urna pretium. Donec luctus non leo "
                    "sed euismod. Phasellus in diam et leo fringilla "
                    "adipiscing ullamcorper nec sapien. Sed condimentum metus "
                    "at lacus vehicula vulputate. Nam fermentum faucibus ipsum "
                    "ut gravida. In sed felis tellus. Aliquam imperdiet, augue "
                    "et eleifend cursus, elit risus accumsan justo, eu aliquam "
                    "quam massa id risus. Donec sagittis orci lorem, a "
                    "vulputate lacus sodales ut. Proin massa massa, aliquet "
                    "vitae felis et, porttitor ornare enim.")
            print >>i_file, "\n".join(textwrap.wrap(text, 80)) + "\n\n"

            # Writing the methods sections
            print >>i_file, template.section("Methods")
            text = ("Plink was used for the data cleanup procedure. The "
                    r"automated script \texttt{pyGenClean} "
                    "version " + prog_version + r"~\cite{pyGenClean} was "
                    r"used to launch the analysis.\\")
            print >>i_file, "\n".join(textwrap.wrap(text, 80)) + "\n\n"

            print >>i_file, "Briefly, the clean up procedure was as follow:"
            print >>i_file, r"\begin{enumerate}" + "\n"
            for step, desc in zip(steps, descriptions):
                step = step.replace("_", r"\_")
                to_print = template.item(desc)
                to_print += " ({})".format(template.texttt(step))
                print >>i_file, "\n".join(textwrap.wrap(to_print)) + "\n"
            print >>i_file, r"\end{enumerate}" + "\n"

            # Including the summaries in the results section
            print >>i_file, template.section("Results")
            for name in summaries:
                full_path = os.path.abspath(name)
                if os.path.isfile(full_path):
                    print >>i_file, r"\input{" + full_path + "}\n"

            print >>i_file, r"\begin{thebibliography}{9}"
            entry = template.bib_entry(name="pyGenClean",
                                       authors=("Lemieux Perreault LP, "
                                                "Provost S, Legault MA, "
                                                r"Barhdadi A, Dub\'e MP"),
                                       title=("pyGenClean: efficient tool "
                                              "for genetic data clean up "
                                              "before association testing"),
                                       journal="Bioinformatics",
                                       year="2013",
                                       volume="29",
                                       number="13",
                                       pages="1704--1705")
            print >>i_file, entry
            print >>i_file, r"\end{thebibliography}" + "\n"
            print >>i_file, r"\end{document}"

    except IOError:
        msg = "{}: could not create report".format(report_name)
        raise ProgramError(msg)
