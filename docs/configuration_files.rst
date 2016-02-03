.. _config_files:

Configuration Files
*******************

Two default configuration files are available to run the proposed protocol.
Before using them, be sure to follow the :ref:`preprocessing_label` described in
the :ref:`proposed_protocol_label` section.

Note that lines starting with a ``#`` are comments, and are not used by the
pipeline. The default parameters were commented out, and could be uncommented to
change their values.


.. _first_conf_file:

First Configuration File
========================

This file should be use with the original dataset as input. Only change the
``loop-assoc`` file name in the ``plate_bias`` section (``[8]``) and the
reference population files (``ceu-bfile``, ``yri-bfile`` and ``jpt-chb-bfile``
in the ``check_ethnicity`` section (``[10]``). Those last three datasets are
provided and can be downloaded at `http://www.statgen.org
<http://www.statgen.org>`_.

If you want to generate the gender and BAF and LRR plots, you will require to
provide the intensities (``sex-chr-intensities``  and ``lrr-baf-raw-dir`` in the
``sex_check`` section (``[7]``) after uncommenting the required options).

.. literalinclude:: _static/configuration_files/configuration_example_1_of_2.ini
    :linenos:
    :language: lighttpd


.. _second_conf_file:

Second Configuration File
=========================

This configuration file should be run after the :ref:`first_conf_file` and with
the output of the second sample missingness section (``[6]`` in the
:ref:`first_conf_file`).

A file containing the samples and markers to be removed should be created using
the output of the ``sex_check``, ``find_related_samples``, ``check_ethnicity``
and ``plate_bias`` sections of the :ref:`first_conf_file`.

.. literalinclude:: _static/configuration_files/configuration_example_2_of_2.ini
    :linenos:
    :language: lighttpd

