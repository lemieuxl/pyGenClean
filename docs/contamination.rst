.. _contamination_label:

Contamination Module
====================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_check_contamination --help
    usage: pyGenClean_check_contamination [-h] [-v] --bfile FILE --raw-dir DIR
                                          [--colsample COL] [--colmarker COL]
                                          [--colbaf COL] [--colab1 COL]
                                          [--colab2 COL] [--sge]
                                          [--sge-walltime TIME]
                                          [--sge-nodes INT INT]
                                          [--sample-per-run-for-sge INT]
                                          [--out FILE]

    Check BAF and LogR ratio for data contamination.

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit

    Input File:
      --bfile FILE          The input file prefix (will find the plink binary
                            files by appending the prefix to the .bim, .bed and
                            .fam files, respectively).

    Raw Data:
      --raw-dir DIR         Directory containing the raw data (one file per
                            sample, where the name of the file (minus the
                            extension) is the sample identification number.
      --colsample COL       The sample column. [default: Sample Name]
      --colmarker COL       The marker column. [default: SNP Name]
      --colbaf COL          The B allele frequency column. [default: B Allele
                            Freq]
      --colab1 COL          The AB Allele 1 column. [default: Allele1 - AB]
      --colab2 COL          The AB Allele 2 column. [default: Allele2 - AB]

    SGE Options:
      --sge                 Use SGE for parallelization.
      --sge-walltime TIME   The walltime for the job to run on the cluster. Do not
                            use if you are not required to specify a walltime for
                            your jobs on your cluster (e.g. 'qsub
                            -lwalltime=1:0:0' on the cluster).
      --sge-nodes INT INT   The number of nodes and the number of processor per
                            nodes to use (e.g. 'qsub -lnodes=X:ppn=Y' on the
                            cluster, where X is the number of nodes and Y is the
                            number of processor to use. Do not use if you are not
                            required to specify the number of nodes for your jobs
                            on the cluster.
      --sample-per-run-for-sge INT
                            The number of sample to run for a single SGE job.
                            [default: 30]

    Output File:
      --out FILE            The prefix of the output files. [default:
                            contamination]


Input Files
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

For more information about the actual algorithms and source codes, refer to the
following page.

* :py:mod:`pyGenClean.Contamination.contamination`
