.. _proposed_protocol_label:

Proposed Protocol
*****************


.. _contamination_module_label:

Contamination Module
====================

.. note::

    **Input files:**

    * PLINK binary pedfiles (BED, BIM and FAM)

1.  Examine the output file named ``contamination.bafRegress``. It will contain
    the contamination estimates (along with confidence values). Usually, an
    estimate higher than 0.01 means possible contamination.

2.  The automatic report will contain the list of samples with possible
    contamination (*i.e.* an estimate value higher than 0.01).


.. _preprocessing_label:

Preprocessing Steps
===================

*   Remove SNPs without chromosomal and physical position (chromosome and
    position of 0).

*   Remove INDELs (markers with alleles ``I`` or ``D``).

*   Determine if there are duplicated samples. These samples must have exactly
    the same family (``FID``) and individual (``IID``) identification to be
    treated as duplicated samples by the :ref:`dup_sample_module_label`. PLINK's
    option ``--update-ids`` could be used.

*   If input is a transposed pedfile, be sure to use PLINK's option ``--tab`` to
    produce the appropriate file format.

*   For the :ref:`plate_bias`, a text file explaining the plate distribution of
    each sample must be provided using the option ``--loop-assoc`` in the
    configuration file. The following columns are required (in order, without a
    header):

        *   the family identification;
        *   the individual identification;
        *   and the plate identification.

*   Produce parameter files (see the :ref:`config_files` for details about
    parameter file).

*   To launch the analysis consult the section :ref:`how_ro_run`.


.. _dup_sample_module_label:

Duplicated Samples Module
=========================

.. note::

    **Input files:**

    *   PLINK transposed pedfiles from the :ref:`preprocessing_label`.

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``dup_samples.diff`` to evaluate if some samples have many
    discordant genotypes (this could indicate a possible samples mix up). To
    identify discordant samples, use the following command line:

    .. code-block:: console

        $ cut -f4 dup_samples.diff | sort -k1,1 | uniq -c
        $ cut -f5 dup_samples.diff | sort -k1,1 | uniq -c

    If samples are present more than 10,000 times (for 2.5E-6 SNPs) this could
    indicate a sample mix up.

3.  Examine ``dup_samples.not_good_enough`` to determine if samples have a
    concordance rate below the threshold set by the user. These samples **are
    present** in the ``dup_samples.final.tfam`` if they are the chosen ones.

4.  Examine ``dup_samples.summary`` to evaluate completion rate and concordance
    between the replicates of potentially problematic samples.

5.  Examine ``dup_samples.concordance`` file for the problematic samples; this
    could help to determine which sample is the discordant replicate.

6.  If a sample appears problematic rename it and keep it in the analysis to
    determine if it is a duplicate of another sample (mix up) with the related
    sample module.

If necessary, samples present in the ``dup_samples.not_good_enough`` file can be
removed from the data set with the subset module (see the
:ref:`first_subset_label`). If not, proceed to the
:ref:`dup_marker_module_label`).


.. _first_subset_label:

First Subset Module (optional)
==============================

.. note::

    **Input files:**

    From the :ref:`dup_sample_module_label`:

    *   ``dup_samples.final.tfam``
    *   ``dup_samples.final.tped``

1.  Extract the family (``FID``) and individual (``IID``) identification from
    ``dup_samples.not_good_enough`` with the following command line:

    .. code-block:: console

        $ cut -f3,4 dup_samples.not_good_enough | sed "1d" | sort -k1,1 \
        >     | uniq > samples_to_remove


.. _dup_marker_module_label:

Duplicated Markers Module
==========================

.. note::

    **Input files:**

    From the :ref:`dup_sample_module_label`:

    *   ``dup_samples.final.tfam``
    *   ``dup_samples.final.tped``

    or from the :ref:`first_subset_label`:

    *   ``subset.bed``
    *   ``subset.bim``
    *   ``subset.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script

2.  Examine ``dup_snps.duplicated_marker_names`` to detect SNPs with exactly the
    same name but mapping to different chromosomal location. (This file is not
    produce if no duplicated marker names are identified).

3.  Determine the number of duplicated SNPs merged (same allele, same frequency,
    etc). SNPs merged were removed and are listed in the file
    ``dup_snps.removed_duplicates``. Number of lines in this file corresponds to
    number of SNPs merged. SNPs not merged and reasons why (*e.g.*
    ``homo_hetero``, ``diff_frequency``, ``homo_flip``, etc.) are present in
    file ``dup_snps.problems``.

4.  SNPs with concordance rate below the threshold are present in
    ``dup_snps.not_good_enough``. To have the list of those SNPs:

    .. code-block:: console

        $ grep -w concordance dup_snps.not_good_enough | cut -f1 \
        >     > SNP_with_low_concordance_rate

If necessary, use the subset option in the configuration file to remove the low
concordance rate SNPs (see the :ref:`second_subset_label`).


.. _second_subset_label:

Second Subset Module (optional)
===============================

.. note::

    **Input files:**

    From the :ref:`dup_marker_module_label`:

    *   ``dup_snps.final.tfam``
    *   ``dup_snps.final.tped``

*   Extract SNPs with concordance rate below the threshold set by the user with
    the command line

    .. code-block:: console

        $ grep -w concordance dup_snps.not_good_enough | cut -f1 \
        >     > SNP_with_low_concordance_rate


.. _clean_noCall_hetero:

Clean No Call and Only Heterozygous Markers Module
==================================================

.. note::

    **Input files:**

    From the :ref:`dup_marker_module_label`:

    *   ``dup_snps.final.tfam``
    *   ``dup_snps.final.tped``

    or from the :ref:`second_subset_label`:

    *   ``subset.bed``
    *   ``subset.bim``
    *   ``subset.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  SNPs removed because they are failed are listed in
    ``clean_noCall_hetero.allFailed``.

3.  SNPs removed because they are all heterozygous are listed in
    ``clean_noCall_hetero.allHetero``.


.. _sample_missingness_01:

Sample Missingness Module  (mind 0.1)
=====================================

.. note::

    **Input files:**

    From the :ref:`clean_noCall_hetero`:

    *   ``clean_noCall_hetero.tfam``
    *   ``clean_noCall_hetero.tped``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine PLINK's log file to detect any problem at this step.

3.  Individuals removed because they did not pass the completion rate threshold
    are listed in ``clean_mind.irem``.


.. _marker_missingness:

Marker Missingness Module
=========================

.. note::

    **Input files:**

    From the :ref:`sample_missingness_01`:

    *   ``clean_mind.bed``
    *   ``clean_mind.bim``
    *   ``clean_mind.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine PLINK's log file to detect any problem at this step.

3.  SNPs removed because they did not pass the completion rate threshold are
    listed in ``clean_geno.removed_snps``.


.. _sample_missingness_02:

Sample Missingness Module (mind 0.02)
=====================================

.. note::

    **Input files:**

    From the :ref:`marker_missingness`:

    *   ``clean_geno.bed``
    *   ``clean_geno.bim``
    *   ``clean_geno.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine PLINK's log file to detect any problem at this step.

3.  Individuals removed because they did not pass the completion rate threshold
    are listed in ``clean_mind.irem``.


.. _sex_check:

Sex Check Module
================

.. note::

    **Input files:**

    From :ref:`sample_missingness_02`:

    *   ``clean_mind.bed``
    *   ``clean_mind.bim``
    *   ``clean_mind.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine PLINK's log file to detect any problem at this step.

3.  Examine ``sexcheck.list_problem_sex``, it contains all individuals
    identified by PLINK as having gender problem.

4.  Examine ``sexcheck.chr23_recodeA.raw.hetero`` to determine heterozygosity on
    the X chromosome of problematic samples. Consanguineous females may have low
    heterozygosity on the X chromosome. If many genotyped SNPs are rare,
    heterozygosity may also be low.

5.  Examine ``sexcheck.chr24_recodeA.raw.noCall`` to determine the number of Y
    markers with missing calls. Females have low number of genotypes for Y
    chromosome markers (high values of missing calls), but is often not equal to
    0 probably because some Y markers come from pseudo autosomal regions. Column
    ``nbGeno`` is the total number of genotypes check and ``nbNoCall`` is the
    number of genotypes with missing calls on chromosome Y. Males should have
    low values in this column while females have higher number of missing calls
    but not equal to the total number of genotypes tested.


If probe intensities from X and Y chromosomes are available and the gender plot
has been created:

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``sexcheck.png`` to detect any individuals in the XXY or X0 regions,
    females in the male cluster and males in the female cluster (see the
    :ref:`gender_plot_figure` figure). Confirm if possible the gender problems
    identified with the previous sex check problem step.

If intensities file for each sample are available and the BAF and LRR plot has
been created:

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``sexcheck_sample-id_lrr_baf.png`` for each sample. Usually, females
    have LRR values around 0 (between -0.5 and 0.5) while males have LRR values
    between -0.5 and -1. Females have three lines on BAF graphics: one at 1
    (homozygous for the B allele), one at 0.5 (heterozygous AB) and one at 0
    (homozygous for the A allele). Males have two lines: one at 1 (homozygous
    for the B allele) and one a at 0 (homozygous for the A allele). For more
    details, see the :ref:`sex_check_plots` section of the
    :ref:`sexcheck_module_lable`.

Keep individuals identified with gender problem until the :ref:`related_samples`
(mix up of samples could be resolved at this step).


.. _plate_bias:

Plate Bias Module
=================

.. note::

    **Input files:**

    From the :ref:`sample_missingness_02`:

    *   ``clean_mind.bed``
    *   ``clean_mind.bim``
    *   ``clean_mind.fam``

    or if subset option is used to remove SNPs from ``nof`` file (see below):

    *   ``subset.bed``
    *   ``subset.bim``
    *   ``subset.fam``

1.  Verify if there is a ``nof`` file produce by PLINK when the input files for
    this step were produced (from the :ref:`sample_missingness_02`). The ``nof``
    contains SNPs with no founder genotypes observed. If so, remove the SNPs
    present in the ``nof`` file using the subset tool before launching the plate
    bias analysis.  Those SNPs, if they are not removed will produced an error
    message when PLINK performs the ``loop-assoc`` analysis and the following
    message will be present in PLINK's log file ``plate_bias.log``: "``ERROR:
    FEXACT error 3``". SNPs on chromosome 24 could also produce this error.

2.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

3.  Examine ``plate_bias.log`` to detect any problem at this step.

4.  The ``plate_bias.significant_SNPs.txt`` file contains a list of SNPs with P
    value below the threshold. Care should be taken with those SNPs if
    significant results are obtained in association tests. These SNPs are NOT
    removed from the data set, only flagged.

5.  Low MAF can explain part of plate bias. Examine the output file
    ``plate_bias.significant_SNPs.frq`` to determine if SNPs have low MAF. Other
    reasons explaining plate bias are relatedness or ethnicity of individuals
    assign to the same plates and none of them on other plates.


.. _related_samples:

Related Samples Module
======================

.. note::

    **Input files:**

    From the :ref:`sample_missingness_02`:

    *   ``clean_mind.bed``
    *   ``clean_mind.bim``
    *   ``clean_mind.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  File ``ibs.pruning_0.1.prune.in`` contains the list of uncorrelated SNPs
    used for the IBS analysis

3.  Examine ``ibs.related_individuals_z1.png`` and
    ``ibs.related_individuals_z2.png`` to detect if there are samples in the
    parent-child, duplicated samples, first degree relative and second degree
    relative areas. (see :ref:`ibs_z1_figure` and :ref:`ibs_z2_figure` plots).

4.  File ``ibs.related_individuals`` lists pairs of related individuals. Index
    column indicates group of related samples. Status column indicated the
    probable link between pair of individuals based on :math:`Z_0`, :math:`Z_1`
    and :math:`Z_2` values (see the :ref:`z_values_table` table [for which
    :math:`Z` values are approximation] or
    :py:func:`RelatedSamples.find_related_samples.extractRelatedIndividuals`
    function for thresholds).

5.  If there are known duplicated samples, examine ``ibs.related_individuals``
    to determine if they were identified correctly, if not this could indicate a
    possible samples mix up.

6.  File ``ibs.choosen_related_individuals`` contains a list of related samples
    to keep. One related sample from the pair is randomly selected. If there are
    a group of related individuals, one sample in randomly selected from the
    group. All non selected samples are listed in
    ``ibs.discarded_related_individuals`` and should be removed from the
    analysis at a later stage.

.. _z_values_table:
.. table:: IBD allele sharing values

    +--------------+--------------+-------------+--------------+------------------------------------+
    | Relationship | :math:`k_0`  | :math:`k_1` | :math:`k_2`  | Coancestry                         |
    |              |              |             |              | :math:`\theta = 1/2 k_2 + 1/4 k_1` |
    +==============+==============+=============+==============+====================================+
    | Unrelated    | 1            | 0           | 0            | 0                                  |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Identical    | 0            | 0           | 1            | :math:`1/2`                        |
    | twins        |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Parent-child | 0            | 1           | 0            | :math:`1/4`                        |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Full         | :math:`1/4`  | :math:`1/2` | :math:`1/4`  | :math:`1/4`                        |
    | siblings     |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Half         | :math:`1/2`  | :math:`1/2` | 0            | :math:`1/8`                        |
    | siblings     |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Uncle        | :math:`1/2`  | :math:`1/2` | 0            | :math:`1/8`                        |
    | nephew       |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Grandparent  | :math:`1/2`  | :math:`1/2` | 0            | :math:`1/8`                        |
    | grandchild   |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | Double first | :math:`9/16` | :math:`3/8` | :math:`1/16` | :math:`1/8`                        |
    | cousins      |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+
    | First        | :math:`3/4`  | :math:`1/4` | 0            | :math:`1/16`                       |
    | cousins      |              |             |              |                                    |
    +--------------+--------------+-------------+--------------+------------------------------------+


.. _ethnicity:

Ethnicity Module
================

.. note::

    **Input files:**

    From the :ref:`sample_missingness_02`:

    *   ``clean_mind.bed``
    *   ``clean_mind.bim``
    *   ``clean_mind.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  File ``ethnic.ibs.pruning_0.1.prune.in`` contains the list of uncorrelated
    SNPs used for the MDS analysis.

3.  File ``ethnic.mds.mds`` contains the list of principale components as
    calculated by PLINK.

4.  Examine ``ethnicity.mds.png``, ``ethnicity.before.png``,
    ``ethnicity.after.png`` and ``ethnicity.outliers.png`` to detect samples
    outside the selected cluster (see :ref:`ethnicity_plots_label` generated
    from the :ref:`ethnicity_module_label` for more information).

    If there are too many outliers still present in the data set (*i.e.* radius
    is too large), analysis can be redone using the ``pyGenClean_find_outliers``
    standalone script, using a different value for ``--multiplier``. For more
    information, refer to the :ref:`ethnicity_find_outliers` section of the
    :ref:`ethnicity_module_label`.

5.  Samples outside the selected cluster are listed in ``ethnicity.outliers``.
    If necessary those samples could be removed at a later stage with the subset
    option.


.. _third_subset_label:

Third Subset Module
===================

.. note::

    **Input files:**

    From the :ref:`sample_missingness_02`:

    *   ``clean_mind.bed``
    *   ``clean_mind.bim``
    *   ``clean_mind.fam``

Use the subset module to remove samples with gender problems (the
:ref:`sex_check`), outliers from the ethnicity cluster (the :ref:`ethnicity`),
related samples (the :ref:`related_samples`) and any other samples that need to
be removed from the data set.

*   To produces a file containing all the samples to remove from the dataset:

    .. code-block:: console

        $ cat sexcheck.list_problem_sex_ids ibs.discarded_related_individuals \
        >     ethnicity.outliers > samples_to_remove.txt

    One sample may be removed for more than one reason, hence be present more
    than one time in the final ``samples_to_remove.txt`` file. This is not an
    issue for this step.


.. _heterozygote_haploid:

Heterozygote Haploid Module
===========================

.. note::

    **Input files:**

    From the :ref:`third_subset_label`:

    *   ``subset.bed``
    *   ``subset.bim``
    *   ``subset.fam``

Samples with gender problems **must have been removed before** performing this
module.

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``without_hh_genotypes.log`` to detect any problem at this step.

Number of heterozygous haploid genotypes set to missing are indicated in
``without_hh_genotypes.log`` file.


.. _maf:

Minor Allele Frequency of Zero Module
=====================================

.. note::

    **Input files:**

    From the :ref:`heterozygote_haploid`:

    *   ``without_hh_genotypes.bed``
    *   ``without_hh_genotypes.bim``
    *   ``without_hh_genotypes.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``flag_maf_0.log`` to detect any problem at this step.

3.  The file ``flag_maf_0.na_list`` contains a list of SNPs with minor allele
    frequency of 0.

If necessary, use subset module to remove SNPs with minor allele frequency of 0,
since they were only flagged using the :ref:`fourth_subset_label`.


.. _fourth_subset_label:

Fourth Subset Module (optional)
===============================

.. note::

    **Input files:**

    From the :ref:`heterozygote_haploid`:

    *   ``without_hh_genotypes.bed``
    *   ``without_hh_genotypes.bim``
    *   ``without_hh_genotypes.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``subset.log`` to detect any problem at this step.


.. _hwe:

Hardy Weinberg Equilibrium Module
=================================

.. note::

    **Input files:**

    From the :ref:`heterozygote_haploid`:

    *   ``without_hh_genotypes.bed``
    *   ``without_hh_genotypes.bim``
    *   ``without_hh_genotypes.fam``

    or from the :ref:`fourth_subset_label`:

    *   ``subset.bed``
    *   ``subset.bim``
    *   ``subset.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``flag_hw.threshold_1e-4.log`` and
    ``flag_hw.threshold_Bonferroni.log`` to detect any problem at this step.

3.  The files ``flag_hw.snp_flag_threshold_Bonferroni`` and
    ``flag_hw.snp_flag_threshold_1e-4`` contain  lists of SNPs with P value
    below Bonferroni and below :math:`1 \times 10^{-4}` threshold, respectively.

The markers are only flagged using this module. If you want to remove those
markers, have a look at the :ref:`fifth_subset_label`.


.. _fifth_subset_label:

Fifth Subset Module (optional)
==============================

.. note::

    **Input files:**

    From the :ref:`heterozygote_haploid`:

    *   ``without_hh_genotypes.bed``
    *   ``without_hh_genotypes.bim``
    *   ``without_hh_genotypes.fam``

    or from the :ref:`fourth_subset_label`:

    *   ``subset.bed``
    *   ``subset.bim``
    *   ``subsert.fam``

1.  Examine the log to confirm options used and to detect any problems occurring
    while running the script.

2.  Examine ``subset.log`` to detect any problem at this step.
