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

from string import Template


# The document preamble
preamble = Template(r"""\documentclass[10pt,twoside,english]{scrartcl}

\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}

\usepackage[letterpaper]{geometry}
\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm,footskip=1cm}

\usepackage{fancyhdr}
\pagestyle{fancy}

\setlength{\parindent}{0cm}

\usepackage{color}

\usepackage{babel}

\usepackage[unicode=true,
 bookmarks=true,bookmarksnumbered=false,bookmarksopen=false,
 breaklinks=false,pdfborder={0 0 0},backref=false,colorlinks=true]
 {hyperref}
\hypersetup{pdftitle={Automatic Sequencing Report},
 pdfauthor={Louis-Philippe Lemieux Perreault},
 linkcolor={link_blue},citecolor={link_blue},urlcolor={link_blue}}

\makeatletter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User specified LaTeX commands.
\usepackage[activate={true,nocompatibility},final,tracking=true,kerning=true,spacing=true,factor=1100,stretch=10,shrink=10]{microtype}
% activate={true,nocompatibility} - activate protrusion and expansion
% final - enable microtype; use "draft" to disable
% tracking=true, kerning=true, spacing=true - activate these techniques
% factor=1100 - add 10% to the protrusion amount (default is 1000)
% stretch=10, shrink=10 - reduce stretchability/shrinkability (default is 20/20)

% Some color definitions
\usepackage{xcolor}
\definecolor{blue}{HTML}{0099CC}
\definecolor{red}{HTML}{CC0000}
\definecolor{green}{HTML}{669900}
\definecolor{link_blue}{HTML}{21759B}

% Caption and subcaptions
\usepackage[bf]{caption}
\usepackage{subcaption}

% We want images
\usepackage{float}
\usepackage{graphicx}

% Enables page in landscape mode
\usepackage{pdflscape}

% Color in tables + dashed line
\usepackage{colortbl}
\usepackage{arydshln}

% The Report number...
\newcommand{\ProjectName}{${project_name}}

% Fancy header
\pagestyle{fancy}
\fancyhf{}
\fancyfoot[RO,LE]{{\footnotesize Automatic QC Report}}
\fancyfoot[LO,RE]{{\footnotesize\ProjectName}}
\fancyfoot[CO,CE]{{\footnotesize\thepage}}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{1pt}

% The plain fancy header
\fancypagestyle{plain}{
    \fancyhf{} % clear all header and footer fields
    \renewcommand{\headrulewidth}{0pt}
    \renewcommand{\footrulewidth}{0pt}
}

% Table's numerotations are in Roman numbers
\renewcommand{\thetable}{\Roman{table}}

% Subfigure numbering
\renewcommand{\thesubfigure}{\Alph{subfigure}}

\let\tempone\itemize
\let\temptwo\enditemize
\renewenvironment{itemize}{\tempone\setlength{\itemsep}{0pt}}{\temptwo}

\let\ttempone\enumerate
\let\ttemptwo\endenumerate
\renewenvironment{enumerate}{\ttempone\setlength{\itemsep}{0pt}}{\ttemptwo}

\let\tttempone\description
\let\tttemptwo\enddescription
\renewenvironment{description}{\tttempone\setlength{\itemsep}{0pt}}{\tttemptwo}

\makeatother

\begin{document}

% The title page
\title{Genetic Data Clean Up}
\subtitle{\ProjectName}
\author{\includegraphics[width=0.3\textwidth]{${logo_path}}}
\date{\today}
\maketitle

% The table of contents
\microtypesetup{protrusion=false} % disables protrusion locally in the document
\tableofcontents{}
\listoffigures{}
\microtypesetup{protrusion=true} % enables protrusion
\cleardoublepage
\newpage
""")


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
