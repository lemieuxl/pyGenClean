.. pyGenClean documentation master file, created by
   sphinx-quickstart on Wed Dec 12 15:08:13 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. raw:: latex

    \listoffigures
    \listoftables

Welcome to pyGenClean's documentation!
**************************************

Introduction
============

Genetic association studies making use of high throughput genotyping arrays
need to process large amounts of data in the order of millions of markers per
experiment. The first step of any analysis with genotyping arrays is typically
the conduct of a thorough data clean up and quality control to remove poor
quality genotypes and generate metrics to inform and select individuals for
downstream statistical analysis.

:py:mod:`pyGenClean` is an informatics tool to facilitate and standardize the
genetic data clean up pipeline with genotyping array data. In conjuction with a
source batch-queuing system, the tool minimizes data manipulation errors, it
accelerates the completion of the data clean up process and it provides
informative graphics and metrics to guide decision making for statistical
analysis.

:py:mod:`pyGenClean` is a command tool working on both Linux and Windows
operating systems. Its usage is shown below:

.. code-block:: console

   $ run_pyGenClean --help
   usage: run_pyGenClean [-h] [-v] [--bfile FILE] [--tfile FILE] [--file FILE]
                         [--report-author AUTHOR] [--report-number NUMBER]
                         [--report-background BACKGROUND] --conf FILE

   Runs the data clean up (version 1.7.0).

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
     --report-author AUTHOR
                           The current project number. [default: pyGenClean]
     --report-number NUMBER
                           The current project author. [default: Simple Project]
     --report-background BACKGROUND
                           Text of file containing the background section of the
                           report.

   Configuration File:
     --conf FILE           The parameter file for the data clean up.


The tool consists of multiple standalone scripts that are linked together via a
main script (``run_pyGenClean``) and a configuration file (the ``--conf``
option), the latter facilitating user customization.

The :ref:`proposed_protocol_figure` shows the proposed data cleanup pipeline.
Each box represents a customizable standalone script with a quick description
of its function. Optional manual checks for go-no-go decisions are indicated.

.. _proposed_protocol_figure:

.. figure:: _static/images/Data_clean_up_schema.png
    :align: center
    :width: 100%
    :alt: Data clean up schema.

    Data clean up protocol schema


Installation
============

:py:mod:`pyGenClean` is a Python package that works on both Linux and Windows
operating systems. It requires a set of Python dependencies and PLINK. Complete
installation procedures are available for both Linux (32 and 64 bits) and
Windows in the following sections.

.. toctree::
    :maxdepth: 1

    install_linux
    install_windows


Input Files
===========

To use :py:mod:`pyGenClean`, two sets of files are required: a set of genotype
files and a configuration file.

Genotype Files
--------------

The input files of the main program (``run_pyGenClean``) is one of the
following:

*   PLINK's pedfile format (use :py:mod:`pyGenClean`'s ``--file`` option)
    consist of two files with the following extensions: ``PED`` and ``MAP``.

*   PLINK's transposed pedfile format (use :py:mod:`pyGenClean`'s ``--tfile``
    option) consist of two files with the following extensions: ``TPED`` and
    ``TFAM``.

*   PLINK's binary pedfile format (use :py:mod:`pyGenClean`'s ``--bfile``
    option) consist of three files with the following extensions: ``BED``,
    ``BIM`` and ``FAM``.

For more information about these file formats, have a look at PLINK's website,
in the *Basic usage/data formats* section
(`http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
<http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml>`_).

.. warning::

    If the format used is the *transposed* one, the columns **must** be
    separated using **tabulations**, but alleles of each markers need to be
    separated by a single space.

    To create this exact transposed pedfile format, you need to use the
    following PLINK's options:

    *   ``--recode`` to recode the file.
    *   ``--transposed`` to create an output file in the transposed pedfile
        format.
    *   ``--tab`` to use tabulations.

Configuration File
------------------

To customized :py:mod:`pyGenClean`, a basic configuration file is required. It
tells which script to use in a specific order. It also sets the different
options and input files, so that the analysis is easy to replicate or modify.

The configuration file consists of sections, led by a ``[section]`` header
(contiguous integers which gives the order of the pipeline) and followed by
customization of this particular part of the pipeline. Lines preceded by a
``#`` are comments and are not read by :py:mod:`pyGenClean`.

The following example first removes samples with a missing rate of 10% and
more, then removes markers with a missing rate of 2% and more. Finally, it
removes the samples with a missing rate of 2% and more.

.. code-block:: lighttpd
    :linenos:

    [1]
    # Removes sample with a missing rate higher than 10%.
    script = sample_missingness
    mind = 0.1

    [2]
    # Removes markers with a missing rate higher than 2%.
    script = snp_missingness
    geno = 0.02

    [3]
    # Removes sample with a missing rate higher than 2%.
    script = sample_missingness
    mind = 0.02

For a more thorough example, complete configuration files are available for
download at `http://www.statgen.org <http://www.statgen.org>`_ and are
explained in the :ref:`config_files` section. For a list of available modules
and standalone script, refer to the :ref:`list_of_scripts`.

.. toctree::
    :maxdepth: 1

    list_of_options


Information about the Protocol
==============================

The following sections describe the proposed protocol and provides information
about which file should be looked for quality control. Finally, configuration
files (with all available parameters) are given.

.. toctree::
    :maxdepth: 2

    protocol
    configuration_files
    how_to_run


The algorithm
=============

All the functions used for this project are shown and explained in the
following section:

.. toctree::
    :maxdepth: 1

    data_clean_up_module
    duplicated_samples
    duplicated_markers
    clean_noCall_hetero_snps
    sample_missingness
    marker_missingness
    sex_check
    plate_bias
    remove_heterozygous_haploid
    find_related_samples
    check_ethnicity
    flag_maf_zero
    flag_hw
    compare_gold_standard
    plink_utils


Appendix
========

.. toctree::
    :maxdepth: 1

    result_table
