Duplicated Samples Module
=========================

The usage of the standalone module is shown below:

.. code-block:: console

    $ pyGenClean_duplicated_samples --help
    usage: pyGenClean_duplicated_samples [-h] --tfile FILE
                                         [--sample-completion-threshold FLOAT]
                                         [--sample-concordance-threshold FLOAT]
                                         [--out FILE]

    Extract and work with duplicated samples.

    optional arguments:
      -h, --help            show this help message and exit

    Input File:
      --tfile FILE          The input file prefix (will find the tped and tfam
                            file by appending the prefix to .tped and .tfam,
                            respectively.) The duplicated samples should have the
                            same identification numbers (both family and
                            individual ids.)

    Options:
      --sample-completion-threshold FLOAT
                            The completion threshold to consider a replicate when
                            choosing the best replicates and for creating the
                            composite samples. [default: 0.9]
      --sample-concordance-threshold FLOAT
                            The concordance threshold to consider a replicate when
                            choosing the best replicates and for creating the
                            composite samples. [default: 0.97]

    Output File:
      --out FILE            The prefix of the output files. [default: dup_samples]

Input Files
-----------

This module uses PLINK's transposed pedfile format (``tped`` and ``tfam``
files). For this step to work, the duplicated samples must have the same
identification (family and sample ID). One should keep a file containing the
original identifications before modifying the dataset accordingly.

Procedure
---------

Here are the steps performed by the module:

1.  Reads the ``tfam`` file to find duplicated samples.
2.  Separates the duplicated samples from the unique samples.
3.  Writes the unique samples into a file.
4.  Reads the ``tped`` file and write the pedigree file for the unique samples.
    Saves in memory the pedigree for the duplicated samples. Updates the indexes
    of the duplicated samples.
5.  If there are no duplicated samples, simply create the final file. Stop here.
6.  Computes the completion (for each of the duplicated samples) and the
    concordance of each sample pairs.
7.  Prints statistics (concordance and completion).
8.  Prints the concordance matrix for each duplicated samples.
9.  Prints the ``tped`` and the ``tfam`` file for the duplicated samples.
10. Chooses the best of each duplicates (to keep and to complete) according to
    completion and concordance.
11. Creates a unique ``tped`` and ``tfam`` from the duplicated samples by
    completing the best chosen one with the other samples.
12. Creates the final dataset.


Output Files
------------

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [*i.e.* dup_samples]):

1.  No output file is created.
2.  No output file is created.
3.  Only one of the two PLINK's transposed pedfiles is created:

    *   ``dup_samples.unique_samples.tfam``: the ``tfam`` file containing only
        the unique samples from the original dataset.

4.  The second of the two PLINK's transposed pedfiles is created (see previous
    step):

    *   ``dup_samples.unique_samples.tped``: the ``tped`` file containing only
        the unique samples from the original dataset.

5.  If there are not duplicated samples, the final PLINK's transposed pedfiles
    are created (if not, continue tu next step):

    *   ``dup_samples.final``: the ``tfam`` and ``tped`` final files.

6.  One result file is created:

    *   ``dup_samples.diff``: a file containing the differences in the genotypes
        for each pair of duplicated samples. Each line contains the following
        information: ``name`` the name of the marker, ``famID`` the family ID,
        ``indID`` the individual ID, ``dupIndex_1`` the index of the first
        duplicated sample in the original dataset (since the identification of
        each duplicated samples are the same), ``dupIndex_2`` the index of the
        second duplicated sample in the original dataset, ``genotype_1`` and
        ``genotype_2``, the genotype of the first and second duplicated samples
        for the current marker, respectively.

7.  One result file is created:

    *   ``dup_samples.summary``: the completion and summarized concordance of
        each duplicated sample pair. The first two columns (``origIndex`` and
        ``dupIndex`` are the indexes of the duplicated sample in the original
        and duplicated transposed pedfiles, respectively.

8.  One result file is created

    *   ``dup_samples.concordance``: the pairwise concordance of each duplicated
        samples.

9.  One set of PLINK's transposed pedfiles:

    *   ``dup_samples.duplicated_samples``: the dataset containing the
        duplicated samples from the original dataset.

10. Two output files are created:

    *   ``dup_samples.chosen_samples.info``: a list of samples that were chosen
        for completion according to their completion and summarized concordance
        with their duplicates. Again, their indexes in the original and
        duplicated transposed pedfiles are saved (the two first columns).
    *   ``dup_samples.excluded_samples.info``: a list of samples that were not
        chosen for completion according to their completion and summarized
        concordance with their duplicates.

11. Multiple output files are created, along with on set of PLINK's transposed
    pedfiles:

    *   ``dup_samples.zeroed_out``: the list of genotypes that were zeroed out
        during completion of the chosen duplicated samples.
    *   ``dup_samples.not_good_enough``: the list of samples that were not good
        enough (according to completion and concordance) to create the composite
        sample (the chosen duplicated samples).

12. Two sets of PLINK's transposed pedfiles are created:

    *   ``dup_samples.chosen_samples``: a transposed pedfiles containing the
        completed chosen samples.
    *   ``dup_samples.final``: the final dataset.

The Algorithm
-------------

For more information about the actual algorithms and source codes (the
:py:mod:`pyGenClean.DupSamples.duplicated_samples` module), refer to the
following sections.

duplicated_samples
..................

.. automodule:: pyGenClean.DupSamples.duplicated_samples
    :members:
