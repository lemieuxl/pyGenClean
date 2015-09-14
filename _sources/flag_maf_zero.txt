Minor Allele Frequency of Zero Module
=====================================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_flag_maf_zero --help
    usage: pyGenClean_flag_maf_zero [-h] --bfile FILE [--out FILE]
    
    Flag SNPs with MAF of 0.

    optional arguments:
      -h, --help    show this help message and exit

    Input File:
      --bfile FILE  The input file prefix (will find the plink binary files by
                    appending the prefix to the .bim, .bed and .fam files,
                    respectively.

    Output File:
      --out FILE    The prefix of the output files. [default: flag_maf_0]

Input FIles
-----------

This module uses PLINK's binary file format (``bed``, ``bim`` and ``fam`` files)
for the source data set (the data of interest).

Procedure
---------

Here are the steps performed by the module:

1.  Computes the frequencies using Plink.
2.  Finds markers with a MAF of zero.

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* ``flag_maf_0``]):

1.  One file and one set of PLINK's result file:

    *   ``flag_maf_0``: the frequency of each marker in the source dataset.
    *   ``flag_maf_0.list``: the list of markers with a minor allele frequency
        of zero.

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`pyGenClean.FlagMAF.flag_maf_zero` module), refer to the following
sections.

pyGenClean.FlagMAF.flag_maf_zero
................................

.. automodule:: pyGenClean.FlagMAF.flag_maf_zero
    :members:
