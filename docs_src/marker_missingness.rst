Marker Missingness Module
=========================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_snp_missingness --help
    usage: pyGenClean_snp_missingness [-h] --bfile FILE [--geno FLOAT]
                                      [--out FILE]

    Computes sample missingness using Plink

    optional arguments:
      -h, --help    show this help message and exit

    Input File:
      --bfile FILE  The input file prefix (will find the plink binary files by
                    appending the prefix to the .bim, .bed and .fam files,
                    respectively.

    Options:
      --geno FLOAT  The missingness threshold (remove SNPs with more than x
                    percent missing genotypes). [Default: 0.020]

    Output File:
      --out FILE    The prefix of the output files. [default: clean_geno]

Input FIles
-----------

This module uses PLINK's binary file format (``bed``, ``bim`` and ``fam`` files)
for the source data set (the data of interest).

Procedure
---------

Here are the steps performed by the module:

1.  Runs Plink to remove markers with a missing rate above a user defined
    threshold.
2.  Finds the markers that were removed (those that have a missing rate above
    the user defined threshold.

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* ``clean_geno``]):

1.  One set of Plink output files:

    *   ``clean_geno.fam``: the dataset with markers having a high missing rate
        removed (according to a user defined threshold).

2.  One custom file:

    *   ``clean_geno.removed_snps``: a list of markers that have a high missing
        rate (above a user defined threshold).

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`pyGenClean.MarkerMissingness.snp_missingness` module), refer to the
following sections.

snp_missingness
...............

.. automodule:: pyGenClean.MarkerMissingness.snp_missingness
    :members:
