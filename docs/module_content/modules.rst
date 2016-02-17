Main pyGenClean pipeline
========================

Here is the usage of the main pipeline.

.. code-block:: console

    $ run_pyGenClean --help
    usage: run_pyGenClean [-h] [-v] [--bfile FILE] [--tfile FILE] [--file FILE]
                          [--report-title TITLE] [--report-author AUTHOR]
                          [--report-number NUMBER]
                          [--report-background BACKGROUND] --conf FILE

    Runs the data clean up (pyGenClean version 1.8.2).

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit

    Input File:
      --bfile FILE          The input file prefix (will find the plink binary
                            files by appending the prefix to the .bim, .bed and
                            .fam files, respectively).
      --tfile FILE          The input file prefix (will find the plink transposed
                            files by appending the prefix to the .tped and .tfam
                            files, respectively).
      --file FILE           The input file prefix (will find the plink files by
                            appending the prefix to the .ped and .fam files).

    Report Options:
      --report-title TITLE  The report title. [default: Genetic Data Clean Up]
      --report-author AUTHOR
                            The current project number. [default: pyGenClean]
      --report-number NUMBER
                            The current project author. [default: Simple Project]
      --report-background BACKGROUND
                            Text of file containing the background section of the
                            report.

    Configuration File:
      --conf FILE           The parameter file for the data clean up.


Algorithm
---------

For more information about the algorithms, have a look below.

.. toctree::
   :maxdepth: 4

   pyGenClean
