.. _automatic_report:

Automatic Report
****************

The :py:mod:`pyGenClean` module creates an automatic report using
:math:`\LaTeX`. Each of the output directory will contain a final report (named
``automatic_report.tex``). The report, once compiled into the PDF format,
includes various quality metrics. To compile the report, perform the following
commands:

.. code-block:: console

    $ pdflatex automatic_report.tex
    $ pdflatex automatic_report.tex
    $ pdflatex automatic_report.tex

.. note::

    Executing the command three times insure that the links in the automatic
    report are properly formed.


Report merging
==============

On a typical data clean up pipeline, multiple directories will be created (one
for each of the parts of the pipeline. A script is provided to merge all those
reports (one per ``data_clean_up.YYYY-MM-DD_HH.MM.SS`` directory) into a single
report. Here is the usage of this script:

.. code-block:: console

    $ pyGenClean_merge_reports --help
    usage: pyGenClean_merge_reports [-h] --qc-dir DIR [DIR ...]
                                    [--report-title TITLE]
                                    [--report-author AUTHOR]
                                    [--report-number NUMBER]
                                    [--report-background BACKGROUND]
                                    [--out-dir FILE]

    Merges automatic reports from other pyGenClean runs.

    optional arguments:
      -h, --help            show this help message and exit

    Input:
      --qc-dir DIR [DIR ...]
                            A list of directory containing pyGenClean runs.

    Report Options:
      --report-title TITLE  The report title. [default: Genetic Data Clean Up]
      --report-author AUTHOR
                            The current project number. [default: pyGenClean]
      --report-number NUMBER
                            The current project author. [default: Simple Project]
      --report-background BACKGROUND
                            Text of file containing the background section of the
                            report.

    Output Directory:
      --out-dir FILE        The name of the directory that will contain the final
                            report. [default: pyGenClean_report]


To execute the report merging procedure, perform the following command:

.. code-block:: console

    $ pyGenClean_merge_reports --qc-dir data_clean_up.*


.. note::

    It is possible to change the title, the author, the number or the
    background of the automatic report by using the ``--report-title``,
    ``--report-author``, ``--report-number`` or ``--report-background``
    options, respectively.
