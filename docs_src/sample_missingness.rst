Sample Missingness Module
=========================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_sample_missingness --help
    usage: pyGenClean_sample_missingness [-h] --ifile FILE [--is-bfile]
                                         [--mind FLOAT] [--out FILE]

    Computes sample missingness using Plink

    optional arguments:
      -h, --help    show this help message and exit

    Input File:
      --ifile FILE  The input file prefix (by default, this input file must be a
                    tfile. If options --is-bfile is used, the input file must be a
                    bfile).

    Options:
      --is-bfile    The input file (--ifile) is a bfile instead of a tfile.
      --mind FLOAT  The missingness threshold (remove samples with more than x
                    percent missing genotypes). [Default: 0.100]

    Output File:
      --out FILE    The prefix of the output files (wich will be a Plink binary
                    file). [default: clean_mind]

Input FIles
-----------

This module uses either PLINK's binary file format (``bed``, ``bim`` and ``fam``
files) or the transposed pedfile format separated by tabulations (``tped`` and
``tfam``) for the source data set (the data of interest).

Procedure
---------

Here are the steps performed by the module:

1.  Uses Plink to remove samples with a high missing rate (above a user defined
    threshold).

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* ``clean_geno``]):

1.  One set of PLINK's output and result files:

    *   ``clean_mind``: the new dataset with samples having a high missing rate
        removed (above a user defined threshold). The file ``clean_mind.irem``
        contains a list of samples that were removed.

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`pyGenClean.SampleMissingness.sample_missingness` module), refer to the
following sections.

sample_missingness
..................

.. automodule:: pyGenClean.SampleMissingness.sample_missingness
    :members:
