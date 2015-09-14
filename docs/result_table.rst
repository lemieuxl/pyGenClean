Result Summary Table
====================

This table summurized  information available from output files produced by
*pyGenClean* during the data clean up procedure. Numbers correspond to number of
lines in output files see :ref:`proposed_protocol_label` for details. Only
removed SNPs and IDs are indicated in the column SNPs and IDs, flagged SNPs or
IDs are present in the :math:`n` column. 

.. tabularcolumns:: p{11.5cm}rrr
.. _result_summary_table:

.. table:: Summary information of the data clean up procedure.

    +---------------------------------------+---------------+-----------+-------+
    | Description                           | :math:`n`     | SNPs      | IDs   |
    +=======================================+===============+===========+=======+
    | Total number of SNPs in file received | 2,379,855     |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Total number of samples               | 494           |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of duplicate samples           | 0             |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of individuals with no         | 0             |           |       |
    | genotype (failed)                     |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPs with no physical       | 7,239         | -7,239    |       |
    | position (chromosome and physical     |               |           |       |
    | position = 0)                         |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of INDEL                       | 43            | -43       |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of replicate controls          | 5             |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of replicate samples           | 0             |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of duplicate SNPs (by          | 5,643         |           |       |
    | chromosome and physical position)     |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Duplicated SNPs by chromosome and     | 5,417         | -5,147    |       |
    | physical position with the same       |               |           |       |
    | allele (merge)                        |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of duplicated SNP with <98%    | 22            | -22       |       |
    | concordance                           |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Completely failed SNPs                | 1             | -1        |       |
    +---------------------------------------+---------------+-----------+-------+
    | All heterozogous SNPs                 | 0             |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of individuals removed because | 5             |           | -5    |
    | they have more than 10% missing       |               |           |       |
    | genotypes                             |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPs removed because they   | 128,562       | -128,562  |       |
    | have more than 2% missing value       |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of individuals removed because | 7             |           | -7    |
    | they have more than 2% missing        |               |           |       |
    | genotypes                             |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of individuals with gender     | 1             |           |       |
    | problem                               |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPs with plate bias test P | 19            |           |       |
    | value below threshold of              |               |           |       |
    | :math:`1\times10^{-7}`                |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPs used for IBS analysis  | 73,651        |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of duplicates pairs or twin    | 1             |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of related pairs               | 2             |           |       |
    | (including twins)                     |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPs used for MDS analysis  | 80,262        |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of individuals with ethnicity  | 20            |           |       |
    | other than Caucasian as detected by   |               |           |       |
    | MDS analysis                          |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of gender problems             | 1             |           | -1    |
    +---------------------------------------+---------------+-----------+-------+
    | Number of related pairs               | 2             |           | -2    |
    +---------------------------------------+---------------+-----------+-------+
    | Number of caucasian outliers          | 20            |           | -20   |
    +---------------------------------------+---------------+-----------+-------+
    | Number of controls                    | 5             |           | -5    |
    +---------------------------------------+---------------+-----------+-------+
    | Number of heterozygote haploid        | 277,206       |           |       |
    | genotypes set to missing (after       |               |           |       |
    | correction of gender problems)        |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPs with MAF=0             | 602,480       | -602,480  |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPS with HWE test P Value  | 603           |           |       |
    | below threshold of                    |               |           |       |
    | :math:`1\times10^{-4}` and higher     |               |           |       |
    | than Bonferroni threshold             |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | Number of SNPS with HWE test below    | 162           | -162      |       |
    | Bonferroni threshold                  |               |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | **Total number of SNPs**              | **1,635,931** |           |       |
    +---------------------------------------+---------------+-----------+-------+
    | **Total number of samples**           | **454**       |           |       |
    +---------------------------------------+---------------+-----------+-------+

