.. _plate_bias_label:

Plate Bias Module
=================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_plate_bias --help
    usage: pyGenClean_plate_bias [-h] [-v] --bfile FILE --loop-assoc FILE
                                 [--pfilter FLOAT] [--out FILE]

    Checks for plate bias.

    optional arguments:
      -h, --help         show this help message and exit
      -v, --version      show program's version number and exit

    Input File:
      --bfile FILE       The input file prefix (will find the plink binary files
                         by appending the prefix to the .bim, .bed and .fam files,
                         respectively.
      --loop-assoc FILE  The file containing the plate organization of each
                         samples. Must contains three column (with no header):
                         famID, indID and plateName.

    Options:
      --pfilter FLOAT    The significance threshold used for the plate effect.
                         [default: 1.0e-07]

    Output File:
      --out FILE         The prefix of the output files. [default: plate_bias]


Input Files
-----------

This module uses PLINK's binary file format (``bed``, ``bim`` and ``fam`` files)
for the source data set (the data of interest).

Procedure
---------

Here are the steps performed by the module:

1.  Runs the plate bias analysis using Plink.
2.  Extracts the list of significant markers after plate bias analysis.
3.  Computes the frequency of all significant markers after plate bias analysis.

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* ``plate_bias``]):

1.  One set of PLINK's result files:

    *   ``plate_bias``: this includes file like
        ``plate_bias.PLATE_NAME.assoc.fisher`` containing a list of markers that
        were significant after the plate bias analysis.

2.  One file:

    *   ``plate_bias.significant_SNPs.txt``: a file containing a list of markers
        that were significant after the plate bias analysis (contained in all
        the ``*.fisher`` files).

3.  One set of PLINK's result files:

    *   ``plate_bias.significant_SNPs``: the result files containing the
        frequency of all the markers that were significant after the plate bias
        analysis (contained in all the ``*.fisher`` files).

The Algorithm
-------------

For more information about the actual algorithms and source codes, refer to the
following page.

* :py:mod:`pyGenClean.PlateBias.plate_bias`
