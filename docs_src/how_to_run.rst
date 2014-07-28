.. _how_ro_run:

How to Run the Pipeline
=======================

.. warning::

    If you are working in a Linux environment, before doing anything, be sure to
    activate the Python virtual environment (refer to the
    :ref:`python_virtual_environment_label` installation section for more
    information).

    If you are working in a Windows environment, you will need to modify the
    command so that it loads correctly. For example, the following command

    .. code-block:: console

        $ run_pyGenClean \
        >     --conf conf_1.txt \
        >     --tfile /PATH/TO/ORIGINAL/DATASET_PREFIX

    will become

    .. code-block:: none

       > python C:\Python27\Scripts\run_pyGenClean ^
           --conf conf_1.txt ^
           --bfile \PATH\TO\ORIGINAL\DATASET_PREFIX

Modify the :ref:`first_conf_file` so that it suits your needs. After following
the :ref:`preprocessing_label` described in the :ref:`proposed_protocol_label`
section, run the following command:

.. code-block:: console

    $ run_pyGenClean \
    >     --conf conf_1.txt \
    >     --tfile /PATH/TO/ORIGINAL/DATASET_PREFIX

While the protocol is running, check the outputs according to the
:ref:`proposed_protocol_label`. If there are any problems, interrupt the
analysis and make the required modifications. The completed steps can be skipped
by commenting them out, while using the last output dataset as the input one for
the steps that need to be done again.

Once everything was checked, locate the samples and the markers that need to be
removed. For example, if the output directory from the first dataset is
``data_clean_up.YYYY-MM-DD_HH.MM.SS``, the following command will help you:

.. code-block:: console

    $ output_dir=data_clean_up.YYYY-MM-DD_HH.MM.SS
    $ cat $output_dir/7_sex_check/sexcheck.list_problem_sex_ids \
    >     $output_dir/9_find_related_samples/ibs.discarded_related_individuals \
    >     $output_dir/10_check_ethnicity/ethnicity.outliers \
    >     > samples_to_remove.txt

Then, modify the first ``subset`` section in the :ref:`second_conf_file` so that
it reads:

.. code-block:: lighttpd
    :linenos:

    [11]
    script = subset
    remove = samples_to_remove.txt
    exclude = data_clean_up.YYYY-MM-DD_HH.MM.SS/8_plate_bias/plate_bias.significant_SNPs.txt

Once everything was checked, run the following command to finish the data clean
up pipeline:

.. code-block:: console

    $ output_dir=data_clean_up.YYYY-MM-DD_HH.MM.SS
    $ run_pyGenClean \
    >     --conf conf_2.txt \
    >     --bfile $output_dir/6_sample_missingness/clean_mind

If you want to removed the markers that were flagged in the ``flag_maf_zero``
and ``flag_hw`` section, performed the following commands (using the newly
created output directory ``data_clean_up.YYYY-MM-DD_HH.MM.SS``):

.. code-block:: console

    $ output_dir=data_clean_up.YYYY-MM-DD_HH.MM.SS
    $ cat $output_dir/13_flag_maf_zero/flag_maf_0.list \
    >     $output_dir/14_flag_hw/flag_hw.snp_flag_threshold_1e-4 \
    >     > markers_to_exclude.txt
    $ pyGenClean_subset_data \
    >     --ifile $output_dir/14_remove_heterozygous_haploid/without_hh_genotypes \
    >     --is-bfile \
    >     --exclude markers_to_exclude.txt \
    >     --out final_dataset
