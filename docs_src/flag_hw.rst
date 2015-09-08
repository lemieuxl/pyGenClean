Hardy Weinberg Equilibrium Module
=================================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_flag_hw --help
    usage: pyGenClean_flag_hw [-h] --bfile FILE [--hwe FLOAT] [--out FILE]
    
    Flag SNPs with Hardy-Weinberg disequilibrium.

    optional arguments:
      -h, --help    show this help message and exit

    Input File:
      --bfile FILE  The input file prefix (will find the plink binary files by
                    appending the prefix to the .bim, .bed and .fam files,
                    respectively.

    Options:
      --hwe FLOAT   The Hardy-Weinberg equilibrium threshold. [default: 1e-4]

    Output File:
      --out FILE    The prefix of the output files. [default: flag_hw]

Input Files
-----------

This module uses PLINK's binary file format (``bed``, ``bim`` and ``fam`` files)
for the source data set (the data of interest).

Procedure
---------

Here are the steps performed by the module:

1.  Computes the number of markers in the input file.
2.  Computes the Bonferroni threshold.
3.  Runs Plink to find failed markers for HWE with the Bonferroni threshold.
4.  Runs Plink to find failed markers for HWE with the default threshold.
5.  Compares the two marker lists (Bonferroni and default threshold) and finds
    markers that are between the two thresholds.

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* ``flag_hw``]):

1.  No files are created for this step.
2.  No files are created for this step.
3.  One set of PLINK's binary file is created:

    *   ``flag_hw.threshold_X.Xe-X``: the data set containing only the markers
        that pass the HWE test (above the Bonferroni threshold).

4.  One set of PLINK's binary file is created:

    *   ``flag_hw.threshold_1e-4``: the data set containing only the markers
        that pass the HWE test (above the genome wide significance threshold of
        :math:`1 \times 10^{-4}`). This value can be modified at the command
        line.

5.  Three files are created:

    *   ``flag_hw.snp_flag_threshold_X.Xe-X``: the list of markers that failed
        HWE test for the Bonferroni threshold.
    *   ``flag_hw.snp_flag_threshold_1e-4``: the list of markers that failed HWE
        test for the genome wide significance threshold of
        :math:`1 \times 10^{-4}`. This value can be modified at the command line.
    *   ``flag_hw.snp_flag_threshold_between_1e-4-X.Xe-X``: the list of markers
        that failed HWE test at a threshold between the Bonferroni and the
        genome wide significance thresholds, so that you can exclude only the
        ones that have a lower p value than the Bonferroni threshold.

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`FlagHW.flag_hw` module), refer to the following sections.

FlagHW.flag_hw
..............

.. automodule:: FlagHW.flag_hw
    :members:
