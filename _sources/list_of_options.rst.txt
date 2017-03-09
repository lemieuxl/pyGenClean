.. _list_of_scripts:

List of Modules and their Options
*********************************

The following sections show a list the available scripts that can be used in the
configuration file, along with their options for customization.


Contamination
=============

The name to use in the configuration file is ``contamination`` and the
:ref:`contamination_table` table shows its configuration.

.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _contamination_table:

.. table:: List of options for the **contamination** script.

    +------------------------------+------------+-----------------------------+
    | Option                       | Type       | Description                 |
    +==============================+============+=============================+
    | ``--raw-dir``                | ``STRING`` | Directory containing the raw|
    |                              |            | data (one file per sample,  |
    |                              |            | where the name of the file  |
    |                              |            | (minus the extension) is the|
    |                              |            | sample identification       |
    |                              |            | number.                     |
    +------------------------------+------------+-----------------------------+
    | ``--colsample``              | ``STRING`` | The sample column.          |
    +------------------------------+------------+-----------------------------+
    | ``--colmarker``              | ``STRING`` | The marker column.          |
    +------------------------------+------------+-----------------------------+
    | ``--colbaf``                 | ``STRING`` | The B allele frequency      |
    |                              |            | column.                     |
    +------------------------------+------------+-----------------------------+
    | ``--colab1``                 | ``STRING`` | The AB Allele 1 column.     |
    +------------------------------+------------+-----------------------------+
    | ``--colab2``                 | ``STRING`` | The AB Allele 2 column.     |
    +------------------------------+------------+-----------------------------+
    | ``--sge``                    |            | Use SGE for parallelization.|
    +------------------------------+------------+-----------------------------+
    | ``--sge-walltime``           | ``STRING`` | The walltime for the job to |
    |                              |            | run on the cluster. Do not  |
    |                              |            | use if you are not required |
    |                              |            | to specify a walltime for   |
    |                              |            | your jobs on your cluster   |
    |                              |            | (*e.g.*                     |
    |                              |            | '``qsub -lwalltime=1:0:0``' |
    |                              |            | on the cluster).            |
    +------------------------------+------------+-----------------------------+
    | ``--sge-nodes``              | ``INT``    | The number of nodes and the |
    |                              |            | number of processor per     |
    |                              |            | nodes to use (*e.g.*        |
    |                              |            | '``qsub -lnodes=X:ppn=Y``'  |
    |                              |            | on the cluster, where X is  |
    |                              |            | the number of nodes and Y is|
    |                              |            | the number of processor to  |
    |                              |            | use. Do  not use if you are |
    |                              |            | not required to specify the |
    |                              |            | number of nodes  for your   |
    |                              |            | jobs on the cluster.        |
    +------------------------------+------------+-----------------------------+
    | ``--sample-per-run-for-sge`` | ``INT``    | The number of sample to run |
    |                              |            | for a single SGE job.       |
    +------------------------------+------------+-----------------------------+

The name of the standalone script is ``pyGenClean_check_contamination``.


.. _duplicated_samples_options:

Duplicated Samples
==================

The name to use in the configuration file is ``duplicated_samples`` and the
:ref:`duplicated_samples_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _duplicated_samples_table:

.. table:: List of options for the **duplicated_samples** script.

    +------------------------------------+-----------+-------------------------+
    | Option                             | Type      | Description             |
    +====================================+===========+=========================+
    | ``--sample-completion-threshold``  | ``FLOAT`` | The completion          |
    |                                    |           | threshold to consider a |
    |                                    |           | replicate when choosing |
    |                                    |           | the best replicates and |
    |                                    |           | for creating the        |
    |                                    |           | composite samples.      |
    |                                    |           | [default: 0.9]          |
    +------------------------------------+-----------+-------------------------+
    | ``--sample-concordance-threshold`` | ``FLOAT`` | The concordance         |
    |                                    |           | threshold to consider a |
    |                                    |           | replicate when choosing |
    |                                    |           | the best replicates and |
    |                                    |           | for creating the        |
    |                                    |           | composite samples.      |
    |                                    |           | [default: 0.97]         |
    +------------------------------------+-----------+-------------------------+

The name of the standalone script is ``pyGenClean_duplicated_samples``.


Duplicated Markers
==================

The name to use in the configuration file is ``duplicated_snps`` and the
:ref:`duplicated_markers_table` table shows its configuration.

.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _duplicated_markers_table:

.. table:: List of options for the **duplicated_snps** script.

    +---------------------------------+-----------+--------------------------+
    | Option                          | Type      | Description              |
    +=================================+===========+==========================+
    | ``--snp-completion-threshold``  | ``FLOAT`` | The completion threshold |
    |                                 |           | to consider a replicate  |
    |                                 |           | when choosing the best   |
    |                                 |           | replicates and for       |
    |                                 |           | composite creation.      |
    |                                 |           | [default: 0.9]           |
    +---------------------------------+-----------+--------------------------+
    | ``--snp-concordance-threshold`` | ``FLOAT`` | The concordance          |
    |                                 |           | threshold to consider a  |
    |                                 |           | replicate when choosing  |
    |                                 |           | the best replicates and  |
    |                                 |           | for composite creation.  |
    |                                 |           | [default: 0.98]          |
    +---------------------------------+-----------+--------------------------+
    | ``--frequency_difference``      | ``FLOAT`` | The maximum difference   |
    |                                 |           | in frequency between     |
    |                                 |           | duplicated markers       |
    |                                 |           | [default: 0.05]          |
    +---------------------------------+-----------+--------------------------+

The name of the standalone script is ``pyGenClean_duplicated_snps``.


Clean No Call and Only Heterozygous Markers
===========================================

The name to use in the configuration file is ``noCall_hetero_snps`` and there
are no customization possible.

The name of the standalone script is ``pyGenClean_clean_noCall_hetero_snps``.


Sample Missingness
==================

The name to use in the configuration file is ``sample_missingness`` and the
:ref:`sample_missingness_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _sample_missingness_table:

.. table:: List of options for the **sample_missingness** script.

    +------------+-----------+------------------------------------------------+
    | Option     | Type      | Description                                    |
    +============+===========+================================================+
    | ``--mind`` | ``FLOAT`` | The missingness threshold (remove samples with |
    |            |           | more than x percent missing genotypes).        |
    |            |           | [Default: 0.100]                               |
    +------------+-----------+------------------------------------------------+

The name of the standalone script is ``pyGenClean_sample_missingness``.


Marker Missingness
==================

The name to use in the configuration file is ``snp_missingness`` and the
:ref:`snp_missingness_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _snp_missingness_table:

.. table:: List of options for the **snp_missingness** script.

    +------------+-----------+---------------------------------------------+
    | Option     | Type      | Description                                 |
    +============+===========+=============================================+
    | ``--geno`` | ``FLOAT`` | The missingness threshold (remove SNPs with |
    |            |           | more than x percent missing genotypes).     |
    |            |           | [Default: 0.020]                            |
    +------------+-----------+---------------------------------------------+

The name of the standalone script is ``pyGenClean_snp_missingness``.


Sex Check
=========

The name to use in the configuration file is ``sex_check`` and the
:ref:`sex_check_table` table shows its configuration.


.. tabularcolumns:: p{6.3cm}Lp{7.5cm}
.. _sex_check_table:

.. table:: List of options for the **sex_check** script.

    +---------------------------+------------+---------------------------------+
    | Option                    | Type       | Description                     |
    +===========================+============+=================================+
    | ``--femaleF``             | ``FLOAT``  | The female F threshold.         |
    |                           |            | [default: < 0.300000]           |
    +---------------------------+------------+---------------------------------+
    | ``--maleF``               | ``FLOAT``  | The male F threshold.           |
    |                           |            | [default: > 0.700000]           |
    +---------------------------+------------+---------------------------------+
    | ``--nbChr23``             | ``INT``    | The minimum number of markers   |
    |                           |            | on chromosome 23 before         |
    |                           |            | computing Plink's sex check     |
    |                           |            | [default: 50]                   |
    +---------------------------+------------+---------------------------------+
    | ``--gender-plot``         |            | Create the gender plot          |
    |                           |            | (summarized chr Y intensities   |
    |                           |            | in function of summarized chr X |
    |                           |            | intensities) for problematic    |
    |                           |            | samples. Not used by default.   |
    +---------------------------+------------+---------------------------------+
    | ``--sex-chr-intensities`` | ``FILE``   | A file containing alleles       |
    |                           |            | intensities for each of the     |
    |                           |            | markers located on the X and Y  |
    |                           |            | chromosome for the gender plot. |
    +---------------------------+------------+---------------------------------+
    | ``--gender-plot-format``  | ``STRING`` | The output file format for the  |
    |                           |            | gender plot (png, ps, or pdf    |
    |                           |            | formats are available).         |
    |                           |            | [default: png]                  |
    +---------------------------+------------+---------------------------------+
    | ``--lrr-baf``             |            | Create the LRR and BAF plot for |
    |                           |            | problematic samples. Not used   |
    |                           |            | by default.                     |
    +---------------------------+------------+---------------------------------+
    | ``--lrr-baf-raw-dir``     | ``DIR``    | Directory or list of            |
    |                           |            | directories containing          |
    |                           |            | information about every samples |
    |                           |            | (BAF and LRR).                  |
    +---------------------------+------------+---------------------------------+
    | ``--lrr-baf-format``      | ``STRING`` | The output file format for the  |
    |                           |            | LRR and BAF plot (png, ps or    |
    |                           |            | pdf formats are available).     |
    |                           |            | [default: png]                  |
    +---------------------------+------------+---------------------------------+
    | ``--lrr-baf-dpi``         | ``INT``    | The pixel density of the        |
    |                           |            | figure(s) (DPI).                |
    +---------------------------+------------+---------------------------------+

The name of the standalone script is ``pyGenClean_sex_check``. If you want to
redo the BAF and LRR plot or the gender plot, you can use the
``pyGenClean_baf_lrr_plot`` and ``pyGenClean_gender_plot`` scripts,
respectively.


Plate Bias
==========

The name to use in the configuration file is ``plate_bias`` and the
:ref:`plate_bias_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _plate_bias_table:

.. table:: List of options for the **plate_bias** script.

    +------------------+-----------+-----------------------------------------+
    | Option           | Type      | Description                             |
    +==================+===========+=========================================+
    | ``--loop-assoc`` | ``FILE``  | The file containing the plate           |
    |                  |           | organization of each samples. Must      |
    |                  |           | contains three column (with no header): |
    |                  |           | famID, indID and plateName.             |
    +------------------+-----------+-----------------------------------------+
    | ``--pfilter``    | ``FLOAT`` | The significance threshold used for the |
    |                  |           | plate effect. [default: 1.0e-07]        |
    +------------------+-----------+-----------------------------------------+

The name of the standalone script is ``pyGenClean_plate_bias``.


Heterozygous Haploid
====================

The name to use in the configuration file is ``remove_heterozygous_haploid`` and
there are no customization possible.

The name of the standalone script is ``pyGenClean_remove_heterozygous_haploid``.


Related Samples
===============

The name to use in the configuration file is ``find_related_samples`` and the
:ref:`find_related_samples_table` table shows its configuration.


.. tabularcolumns:: p{5.1cm}Lp{7.5cm}
.. _find_related_samples_table:

.. table:: List of options for the **find_related_samples** script.

    +-----------------------------+------------+-------------------------------+
    | Option                      | Type       | Description                   |
    +=============================+============+===============================+
    | ``--genome-only``           |            | Only create the genome file.  |
    |                             |            | Not selected by default.      |
    +-----------------------------+------------+-------------------------------+
    | ``--min-nb-snp``            | ``INT``    | The minimum number of markers |
    |                             |            | needed to compute IBS values. |
    |                             |            | [Default: 10000]              |
    +-----------------------------+------------+-------------------------------+
    | ``--indep-pairwise``        | ``INT``    | Three numbers: window size,   |
    |                             | ``INT``    | window shift and the r2       |
    |                             | ``FLOAT``  | threshold. [default: ['50',   |
    |                             |            | '5', '0.1']]                  |
    +-----------------------------+------------+-------------------------------+
    | ``--maf``                   | ``FLOAT``  | Restrict to SNPs with MAF >=  |
    |                             |            | threshold. [default: 0.05]    |
    +-----------------------------+------------+-------------------------------+
    | ``--ibs2-ratio``            | ``FLOAT``  | The initial IBS2* ratio (the  |
    |                             |            | minimum value to show in the  |
    |                             |            | plot. [default: 0.8]          |
    +-----------------------------+------------+-------------------------------+
    | ``--sge``                   |            | Use SGE for parallelization.  |
    +-----------------------------+------------+-------------------------------+
    | ``--sge-walltime``          | ``STRING`` | The time limit (for clusters).|
    |                             |            | Do not use if you are not     |
    |                             |            | required to specify a walltime|
    |                             |            | for your jobs on your cluster |
    |                             |            | (e.g. ``-lwalltime=1:0:0`` on |
    |                             |            | the cluster). Allow enough    |
    |                             |            | time for proper job           |
    |                             |            | completion.                   |
    +-----------------------------+------------+-------------------------------+
    | ``--sge-nodes``             | ``INT``    | The number of nodes and the   |
    |                             | ``INT``    | number of processor per nodes |
    |                             |            | to use (e.g. ``qsub           |
    |                             |            | -lnodes=X:ppn=Y`` on the      |
    |                             |            | cluster, where X is the number|
    |                             |            | of nodes and Y is the number  |
    |                             |            | of processor to use. Do not   |
    |                             |            | use if you are not required to|
    |                             |            | specify the number of nodes   |
    |                             |            | for your jobs on the cluster. |
    |                             |            | Allow enough ressources for   |
    |                             |            | proper job completion.        |
    +-----------------------------+------------+-------------------------------+
    | ``--line-per-file-for-sge`` | ``INT``    | The number of line per file   |
    |                             |            | for SGE task array.           |
    |                             |            | [default: 100]                |
    +-----------------------------+------------+-------------------------------+

The name of the standalone script is ``pyGenClean_find_related_samples``. Even
though randomly choosing a subset of related samples is done automatically, you
can use the ``pyGenClean_merge_related_samples`` to perform it again.


Ethnicity
=========

The name to use in the configuration file is ``check_ethnicity`` and the
:ref:`check_ethnicity_table` table shows its configuration.


.. tabularcolumns:: p{5.1cm}Lp{7.5cm}
.. _check_ethnicity_table:

.. table:: List of options for the **check_ethnicity** script.

    +-----------------------------+------------+-------------------------------+
    | Option                      | Type       | Description                   |
    +=============================+============+===============================+
    | ``--skip-ref-pops``         |            | Perform the MDS computation,  |
    |                             |            | but skip the three reference  |
    |                             |            | panels.                       |
    +-----------------------------+------------+-------------------------------+
    | ``--ceu-bfile``             | ``FILE``   | The input file prefix (will   |
    |                             |            | find the plink binary files   |
    |                             |            | by appending the prefix to    |
    |                             |            | the .bim, .bed and .fam       |
    |                             |            | files, respectively.) for the |
    |                             |            | CEU population.               |
    +-----------------------------+------------+-------------------------------+
    | ``--yri-bfile``             | ``FILE``   | The input file prefix (will   |
    |                             |            | find the plink binary files   |
    |                             |            | by appending the prefix to    |
    |                             |            | the .bim, .bed and .fam       |
    |                             |            | files, respectively.) for the |
    |                             |            | CEU population.               |
    +-----------------------------+------------+-------------------------------+
    | ``--jpt-chb-bfile``         | ``FILE``   | The input file prefix (will   |
    |                             |            | find the plink binary files   |
    |                             |            | by appending the prefix to    |
    |                             |            | the .bim, .bed and .fam       |
    |                             |            | files, respectively.) for the |
    |                             |            | JPT-CHB population.           |
    +-----------------------------+------------+-------------------------------+
    | ``--min-nb-snp``            | ``FILE``   | The minimum number of markers |
    |                             |            | needed to compute IBS values. |
    |                             |            | [Default: 8000]               |
    +-----------------------------+------------+-------------------------------+
    | ``--indep-pairwise``        | ``INT``    | Three numbers: window size,   |
    |                             | ``INT``    | window shift and the r2       |
    |                             | ``FLOAT``  | threshold. [default: ['50',   |
    |                             |            | '5', '0.1']]                  |
    +-----------------------------+------------+-------------------------------+
    | ``--maf``                   | ``INT``    | Restrict to SNPs with MAF >=  |
    |                             |            | threshold. [default: 0.05]    |
    +-----------------------------+------------+-------------------------------+
    | ``--sge``                   |            | Use SGE for parallelization.  |
    +-----------------------------+------------+-------------------------------+
    | ``--sge-walltime``          | ``STRING`` | The time limit (for clusters).|
    |                             |            | Do not use if you are not     |
    |                             |            | required to specify a walltime|
    |                             |            | for your jobs on your cluster |
    |                             |            | (e.g. ``-lwalltime=1:0:0`` on |
    |                             |            | the cluster). Allow enough    |
    |                             |            | time for proper job           |
    |                             |            | completion.                   |
    +-----------------------------+------------+-------------------------------+
    | ``--sge-nodes``             | ``INT``    | The number of nodes and the   |
    |                             | ``INT``    | number of processor per nodes |
    |                             |            | to use (e.g. ``qsub           |
    |                             |            | -lnodes=X:ppn=Y`` on the      |
    |                             |            | cluster, where X is the number|
    |                             |            | of nodes and Y is the number  |
    |                             |            | of processor to use. Do not   |
    |                             |            | use if you are not required to|
    |                             |            | specify the number of nodes   |
    |                             |            | for your jobs on the cluster. |
    |                             |            | Allow enough ressources for   |
    |                             |            | proper job completion.        |
    +-----------------------------+------------+-------------------------------+
    | ``--ibs-sge-walltime``      | ``STRING`` | The time limit (for clusters) |
    |                             |            | for the IBS jobs. Do not use  |
    |                             |            | if you are not required to    |
    |                             |            | specify a walltime for your   |
    |                             |            | jobs on your cluster (e.g.    |
    |                             |            | ``-lwalltime=1:0:0`` on the   |
    |                             |            | cluster). Allow enough time   |
    |                             |            | for proper job completion.    |
    +-----------------------------+------------+-------------------------------+
    | ``--ibs-sge-nodes``         | ``INT``    | The number of nodes and the   |
    |                             | ``INT``    | number of processor per nodes |
    |                             |            | to use for the IBS jobs (e.g. |
    |                             |            | ``qsub                        |
    |                             |            | -lnodes=X:ppn=Y`` on the      |
    |                             |            | cluster, where X is the number|
    |                             |            | of nodes and Y is the number  |
    |                             |            | of processor to use. Do not   |
    |                             |            | use if you are not required to|
    |                             |            | specify the number of nodes   |
    |                             |            | for your jobs on the cluster. |
    |                             |            | Allow enough ressources for   |
    |                             |            | proper job completion.        |
    +-----------------------------+------------+-------------------------------+
    | ``--line-per-file-for-sge`` | ``INT``    | The number of line per file   |
    |                             |            | for SGE task array.           |
    |                             |            | [default: 100]                |
    +-----------------------------+------------+-------------------------------+
    | ``--nb-components``         | ``INT``    | The number of component to    |
    |                             |            | compute. [default: 10]        |
    +-----------------------------+------------+-------------------------------+
    | ``--outliers-of``           | ``STRING`` | Finds the outliers of this    |
    |                             |            | population. [default: CEU]    |
    +-----------------------------+------------+-------------------------------+
    | ``--multiplier``            | ``FLOAT``  | To find the outliers, we look |
    |                             |            | for more than x times the     |
    |                             |            | cluster standard deviation.   |
    |                             |            | [default: 1.9]                |
    +-----------------------------+------------+-------------------------------+
    | ``--xaxis``                 | ``STRING`` | The component to use for the  |
    |                             |            | X axis. [default: C1]         |
    +-----------------------------+------------+-------------------------------+
    | ``--yaxis``                 | ``STRING`` | The component to use for the  |
    |                             |            | Y axis. [default: C2]         |
    +-----------------------------+------------+-------------------------------+
    | ``--format``                | ``STRING`` | The output file format (png,  |
    |                             |            | ps, pdf, or X11 formats are   |
    |                             |            | available). [default: png]    |
    +-----------------------------+------------+-------------------------------+
    | ``--title``                 | ``STRING`` | The title of the MDS plot.    |
    |                             |            | [default: C2 in function of   |
    |                             |            | C1 - MDS]                     |
    +-----------------------------+------------+-------------------------------+
    | ``--xlabel``                | ``STRING`` | The label of the X axis.      |
    |                             |            | [default: C1]                 |
    +-----------------------------+------------+-------------------------------+
    | ``--ylabel``                | ``STRING`` | The label of the Y axis.      |
    |                             |            | [default: C2]                 |
    +-----------------------------+------------+-------------------------------+
    | ``--create-scree-plot``     |            | Computes Eigenvalues and      |
    |                             |            | creates a scree plot.         |
    +-----------------------------+------------+-------------------------------+
    | ``--scree-plot-title``      | ``STRING`` | The main title of the scree   |
    |                             |            | plot                          |
    +-----------------------------+------------+-------------------------------+

The name of the standalone script is ``pyGenClean_check_ethnicity``. If you want
to redo the outlier detection using a different multiplier, have a look at the
``pyGenClean_find_outliers`` script. If you want to redo any MDS plot, have a
look at the ``pyGenClean_plot_MDS`` script. If you want to compute the
*Eigenvectors* using the ``smartpca`` tool, have a look at the
``pyGenClean_plot_eigenvalues`` script.


Minor Allele Frequency of Zero
==============================

The name to use in the configuration file is ``flag_maf_zero`` and there
are no customization possible.

The name of the standalone script is ``pyGenClean_flag_maf_zero``.


Hardy Weinberg Equilibrium
==========================

The name to use in the configuration file is ``flag_hw`` and the
:ref:`flag_hw_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _flag_hw_table:

.. table:: List of options for the **flag_hw** script.

    +-----------+-----------+-------------------------------------------+
    | Option    | Type      | Description                               |
    +===========+===========+===========================================+
    | ``--hwe`` | ``FLOAT`` | The Hardy-Weinberg equilibrium threshold. |
    |           |           | [default: 1e-4]                           |
    +-----------+-----------+-------------------------------------------+

The name of the standalone script is ``pyGenClean_flag_hw``.


Subsetting the Data
===================

The name to use in the configuration file is ``subset`` and the
:ref:`subset_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _subset_table:

.. table:: List of options for the **subset** script.

    +---------------+----------+--------------------------------------------+
    | Option        | Type     | Description                                |
    +===============+==========+============================================+
    | ``--exclude`` | ``FILE`` | A file containing SNPs to exclude from the |
    |               |          | data set.                                  |
    +---------------+----------+--------------------------------------------+
    | ``--extract`` | ``FILE`` | A file containing SNPs to extract from the |
    |               |          | data set.                                  |
    +---------------+----------+--------------------------------------------+
    | ``--remove``  | ``FILE`` | A file containing samples (FID and IID) to |
    |               |          | remove from the data set.                  |
    +---------------+----------+--------------------------------------------+
    | ``--keep``    | ``FILE`` | A file containing samples (FID and IID) to |
    |               |          | keep from the data set.                    |
    +---------------+----------+--------------------------------------------+

The name of the standalone script is ``pyGenClean_subset_data``.


Comparison with a Gold Standard
===============================

The name to use in the configuration file is ``compare_gold_standard`` and the
:ref:`compare_gold_standard_table` table shows its configuration.


.. tabularcolumns:: p{6.6cm}Lp{7.5cm}
.. _compare_gold_standard_table:

.. table:: List of options for the **compare_gold_standard** script.

    +------------------------+----------+--------------------------------------+
    | Option                 | Type     | Description                          |
    +========================+==========+======================================+
    | ``--gold-bfile``       | ``FILE`` | The input file prefix (will find the |
    |                        |          | plink binary files by appending the  |
    |                        |          | prefix to the .bim, .bed and .fam    |
    |                        |          | files, respectively.) for the Gold   |
    |                        |          | Standard .                           |
    +------------------------+----------+--------------------------------------+
    | ``--same-samples``     | ``FILE`` | A file containing samples which are  |
    |                        |          | present in both the gold standard    |
    |                        |          | and the source panel. One line by    |
    |                        |          | identity and tab separated. For each |
    |                        |          | row, first sample is Gold Standard,  |
    |                        |          | second is source panel.              |
    +------------------------+----------+--------------------------------------+
    | ``--source-manifest``  | ``FILE`` | The illumina marker manifest.        |
    +------------------------+----------+--------------------------------------+
    | ``--source-alleles``   | ``FILE`` | A file containing the source alleles |
    |                        |          | (TOP). Two columns (separated by     |
    |                        |          | tabulation, one with the marker      |
    |                        |          | name, the other with the alleles     |
    |                        |          | (separated by space). No header.     |
    +------------------------+----------+--------------------------------------+
    | ``--sge``              |          | Use SGE for parallelization.         |
    +------------------------+----------+--------------------------------------+
    | ``--do-not-flip``      |          | Do not flip SNPs. WARNING: only use  |
    |                        |          | this option only if the Gold         |
    |                        |          | Standard was generated using the     |
    |                        |          | same chip (hence, flipping is        |
    |                        |          | unnecessary).                        |
    +------------------------+----------+--------------------------------------+
    | ``--use-marker-names`` |          | Use marker names instead of (chr,    |
    |                        |          | position). WARNING: only use this    |
    |                        |          | options only if the Gold Standard    |
    |                        |          | was generated using the same chip    |
    |                        |          | (hence, they have the same marker    |
    |                        |          | names).                              |
    +------------------------+----------+--------------------------------------+

The name of the standalone script is ``pyGenClean_compare_gold_standard``.
