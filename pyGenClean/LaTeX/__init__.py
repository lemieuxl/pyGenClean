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


from string import Template

import jinja2


# The document template
jinja2_env = jinja2.Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%-',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=jinja2.PackageLoader(__name__, "templates"),
)


def section(name):
    """Creates a new section in LaTeX.

    :param name: the name of the new section.
    :type name: string

    :returns: a LaTeX string containing a new section.

    """
    return r"\section{" + name + "}"


def subsection(name):
    """Creates a new section in LaTeX.

    :param name: the name of the new section.
    :type name: string

    :returns: a LaTeX string containing a new section.

    """
    return r"\subsection{" + name + "}"


def item(text):
    """Returns an LaTeX item (for enumerate or itemize).

    :param text: the text.
    :type text: string

    :returns: a LaTeX item.

    """
    return r"\item " + text


def texttt(text):
    """Returns type writer font.

    :param text: the text.
    :type text: string.

    :returns: a type writer font.

    """
    return r"\texttt{" + text + "}"


def textit(text):
    """Returns italic font.

    :param text: the text.
    :type text: string.

    :returns: an italic font.

    """
    return r"\textit{" + text + "}"


def bib_entry(**kwargs):
    """Creates a bibliography entry.

    :param name: the name of the entry.
    :param authors: the authors.
    :param title: the title.
    :param journal: the journal.
    :param year: the year.
    :param volume: the volume.
    :param number: the number.
    :param pages: the pages.

    :type name: string
    :type authors: string
    :type title: string
    :type journal: string
    :type year: string
    :type volume: string
    :type number: string
    :type pages: string

    :returns: a bib entry.

    """
    bib_entry = Template(r"""\bibitem{${name}}
    ${authors}:
    \textbf{${title}.}
    \emph{${journal}} ${year},
    \textbf{${volume}}(${number}): ${pages}""")

    return bib_entry.substitute(kwargs)
