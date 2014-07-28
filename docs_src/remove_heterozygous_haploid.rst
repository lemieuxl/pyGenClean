Heterozygous Haploid Module
===========================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_remove_heterozygous_haploid --help
    usage: pyGenClean_remove_heterozygous_haploid [-h] --bfile FILE [--out FILE]

    Removes heterozygous haploid genotypes.

    optional arguments:
      -h, --help    show this help message and exit

    Input File:
      --bfile FILE  The input file prefix (will find the plink binary files by
                    appending the prefix to the .bim, .bed and .fam files,
                    respectively.

    Output File:
      --out FILE    The prefix of the output files. [default:
                    without_hh_genotypes]

Input Files
-----------

This module uses PLINK's binary file format (``bed``, ``bim`` and ``fam`` files)
for the source data set (the data of interest).

Procedure
---------

Here are the steps performed by the module:

1.  Uses Plink to remove the heterozygous haploid genotypes.

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.*
``without_hh_genotypes``]):

1.  On set of Plink's output files:

    *   ``without_hh_genotypes``: the data set with heterozygous haploid
        genotypes removed.

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`HeteroHap.remove_heterozygous_haploid` module), refer to the following
sections.

HeteroHap.remove_heterozygous_haploid
.....................................

.. automodule:: HeteroHap.remove_heterozygous_haploid
    :members:
