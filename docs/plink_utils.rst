.. _plink_utils_label:

Plink Utils
===========

This module provides useful functions and scripts for efficient interactions
with PLINK's output files. For example, the majority of PLINK's output files are
spaced delimited, and are formated in such a way that it is "beautiful" to the
human eye, but is a bit harder to parse using a script compared to tabulated
files. The :py:func:`pyGenClean.PlinkUtils.createRowFromPlinkSpacedOutput`
function helps producing an array of all the fields for each line.

Comparing BIM files
-------------------

Another example is the fact that when PLINK removes a certain amount of markers
from the data file, it just gives the number of excluded markers, but not a
list. The :py:mod:`pyGenClean.PlinkUtils.compare_bim` module creates a list of
markers that were removed from the original dataset when compared with the new
one. Here is the usage of the standalone script:

.. code-block:: console

    $ pyGenClean_compare_bim --help
    usage: pyGenClean_compare_bim [-h] [-v] --before FILE --after FILE
                                  [--out FILE]

    Compares BIM file.

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit

    Input File:
      --before FILE  The name of the bim FILE before modification.
      --after FILE   The name of the bim FILE after modification.

    Output File:
      --out FILE     The prefix of the output files. [default: snp_removed]


Subsetting a dataset
--------------------

A useful standalone script is the :py:mod:`pyGenClean.PlinkUtils.subset_data`
module. It helps in subsetting a dataset by keeping or removing a set of
samples, and at the same time extracting or excluding a set of markers. The
following standalone script is available for the user:

.. code-block:: console

    $ pyGenClean_subset_data --help
    usage: pyGenClean_subset_data [-h] [-v] --ifile FILE [--is-bfile] [--is-tfile]
                                  [--is-file] [--exclude FILE] [--extract FILE]
                                  [--remove FILE] [--keep FILE] [--out FILE]

    Subsets genotype data using Plink.

    optional arguments:
      -h, --help      show this help message and exit
      -v, --version   show program's version number and exit

    Input File:
      --ifile FILE    The input file prefix. The format will be specified by --is-
                      bfile, --is-tfile or --is-file, for bfile, tfile and file,
                      respectively.
      --is-bfile      The file specified by --ifile is a bfile
      --is-tfile      The file specified by --ifile is a tfile
      --is-file       The file specified by --ifile is a file

    Options:
      --exclude FILE  A file containing SNPs to exclude from the data set.
      --extract FILE  A file containing SNPs to extract from the data set.
      --remove FILE   A file containing samples (FID and IID) to remove from the
                      data set.
      --keep FILE     A file containing samples (FID and IID) to keep from the
                      data set.

    Output File:
      --out FILE      The prefix of the output files. [default: subset]


The standalone script works with the three most used PLINK's format: pedfile,
transposed and binary pedfiles. The ``--is-bfile``, ``--is-tfile`` and
``--is-file`` options tell the standalone script what is the format of the input
file. The output file format will be the same as the input one.

The Algorithm
-------------

For more information about the actual algorithms and source codes, refer to the
following pages.

* :py:mod:`pyGenClean.PlinkUtils`
* :py:mod:`pyGenClean.PlinkUtils.compare_bim`
* :py:mod:`pyGenClean.PlinkUtils.subset_data`
