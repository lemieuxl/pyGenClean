Duplicated Markers Module
=========================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_duplicated_snps --help
    usage: pyGenClean_duplicated_snps [-h] [-v] --tfile FILE
                                      [--snp-completion-threshold FLOAT]
                                      [--snp-concordance-threshold FLOAT]
                                      [--frequency_difference FLOAT] [--out FILE]

    Extracts and merges duplicated markers.

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit

    Input File:
      --tfile FILE          The input file prefix (will find the tped and tfam
                            file by appending the prefix to .tped and .tfam,
                            respectively. A .map file is also required.

    Options:
      --snp-completion-threshold FLOAT
                            The completion threshold to consider a replicate when
                            choosing the best replicates and for composite
                            creation. [default: 0.9]
      --snp-concordance-threshold FLOAT
                            The concordance threshold to consider a replicate when
                            choosing the best replicates and for composite
                            creation. [default: 0.98]
      --frequency_difference FLOAT
                            The maximum difference in frequency between duplicated
                            markers [default: 0.05]

    Output File:
      --out FILE            The prefix of the output files. [default: dup_snps]


Input Files
-----------

This module uses PLINK's transposed pedfile format (``tped`` and ``tfam``
files). It also requires a ``map`` file to speed up the process of finding the
duplicated markers, so that the ``tped`` file is not read.

Procedure
---------

Here are the steps performed by the module:

1.  Reads the ``map`` file to gather marker's position.
2.  Reads the ``tfam`` file.
3.  Finds the unique markers using the ``map`` file.
4.  Process the ``tped`` file, finding unique and duplicated markers according
    to chromosomal positions.
5.  If there are no duplicated markers, stop here.
6.  If there are duplicated markers, print a ``tped`` and ``tfam`` file
    containing the duplicated markers.
7.  Computes the frequency of the duplicated markers (using Plink) and read
    the output file.
8.  Computes the concordance and pairwise completion of each of the
    duplicated markers.
9.  Prints the problematic duplicated markers with a file containing the
    summary of the statistics (completion and pairwise concordance).
10. Print the pairwise concordance in a file (matrices).
11. Choose the best duplicated markers using concordance and completion.
12. Completes the chosen markers with the remaining duplicated markers.
13. Creates the final ``tped`` file, containing the unique markers, the
    chosen duplicated markers that were completed, and the problematic
    duplicated markers (for further analysis). This set excludes markers
    that were used for completing the chosen ones.

Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* dup_snps]):

1.  If the marker names are not unique, one file is created:

    *   ``dup_snps.duplicated_marker_names``: a list of marker names and
        chromosomal positions for each marker with duplicated names. This file
        is not created if there are no markers with duplicated names.

2.  No files are created.

3.  No files are created.

4.  One set of transposed pedfiles.

    *   ``dup_snps.unique_snps``: the transposed pedfiles containing the unique
        markers (according to chromosomal positions).

5.  If there are no duplicated markers (according to chromosomal positions), the
    transposed pedfiles created at the previous step are copied to a new set of
    transposed pedfiles.

    *   ``dup_snps.final``: the final transposed pedfiles.

6.  One set of transposed pedfiles.

    *   ``dup_snps.duplicated_snps``: the transposed pedfiles containing the
        duplicated markers (according to chromosomal positions).

7.  One set of PLINK's result file.

    *   ``dup_snps.duplicated_snps``: the file with the ``frq`` extension
        contains the frequency of each duplicated markers.

8.  No files are created.

9.  Multiple files are created.

    *   ``dup_snps.summary``: contains the completion and pairwise concordance
        between duplicated markers.
    *   ``dup_snps.problems``: contains the list of markers with "problems" that
        can't be used for further completion of the duplicated markers.
        (either a difference in MAF [``diff_frequency``], a difference in the
        minor allele [``diff_minor_allele``], two homozygous markers where one
        is flipped [``homo_flip``], markers with flipped alleles [``flip``], one
        marker is homozygous, the other is heterozygous [``homo_hetero``], one
        marker is homozygous, the other is heterozygous but one is flipped
        [``homo_hetero_flip``] or any other problem [``problem``].

10. One output file is created.

    *   ``dup_snps.concordance``: a matrix containing a pairwise concordance
        comparison for each duplicated markers.

11. Two output files are created.

    *   ``dup_snps.chosen_snps.info``: the list of duplicated markers that were
        chosen for completion with the other markers (the best of the duplicated
        markers, according to concordance and completion).
    *   ``dup_snps.not_chosen_snps.info``: the list of duplicated markers that
        were not chosen for completion with the other markers.

12. Multiple output files are created along with a set of transposed pedfiles.

    *   ``dup_snps.zeroed_out``: the list of genotypes that were zeroed out
        while completing the chosen duplicated markers with the others. Each
        line contains the id of the sample and the name of the marker that was
        zeroed out.
    *   ``dup_snps.not_good_enough``: the list of markers that were not good
        enough (according to concordance and completion) to complete the best of
        the duplicated markers.
    *   ``dup_snps.removed_duplicates``: the list of markers that were used to
        complete the chosen duplicated markers. Those markers were removed from
        the dataset.
    *   ``dup_snps.chosen_snps``: the transposed pedfiles containing the
        completed chosen duplicated markers (a composite of all the duplicated
        markers that were good enough).

13. On set of transposed pedfiles.

    *   ``dup_snps.final``: the final dataset, containing the unique markers,
        the chosen duplicated markers that were complete (composite) and the
        duplicated markers that weren't completed because of various problems.

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`pyGenClean.DupSNPs.duplicated_snps` module), refer to the following
sections.

pyGenClean.DupSNPs.duplicated_snps
..................................

.. automodule:: pyGenClean.DupSNPs.duplicated_snps
    :members:
