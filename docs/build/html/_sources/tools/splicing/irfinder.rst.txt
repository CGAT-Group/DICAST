
.. Links

.. _manual: *not available*
.. |tool| replace:: Irfinder

Irfinder
========

.. note::

  |tool| can be used with bam files or fastq files as reference. Set in :guilabel:`/scripts/asevent_config.sh` the parameter ``use_bam_input_files=1``
  to use bam files and ``use_bam_input_files=0`` to use fastq files.

.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *irfinder*
  **Version:**   *v*
  **License**    *L*
  =============  =================

  **Required files:**

  .. code-block:: bash

    # config.sh
    $gtf
    $fasta
    # use fastq reference
    $fastqdir/*$fastqpair1suffix
    $fastqdir/*$fastqpair2suffix
    #usebamreference
    ???

  **Used parameters**

  .. code-block:: bash

    # config.sh
    $outdir
    # asevent_config.sh
    $use_bam_input_files

|tool|
