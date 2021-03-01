

.. Links

.. _manual: *not available*
.. |tool| replace:: rmats

rmats
=====

.. note::

  |tool| can be used with bam files or fastq files as reference. Set in :guilabel:`/scripts/asevent_config.sh` the parameter ``use_bam_input_files=1``
  to use bam files and ``use_bam_input_files=0`` to use fastq files.


.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *rmats*
  **Version:**   *v*
  **License**    *L*
  =============  =================

  **Required files:**

  .. code-block:: bash

    $gtf
    # use fastq reference
    $fastqdir/*$fastqpair1suffix
    $fastqdir/*$fastqpair2suffix
    #usebamreference
    $controlbam
    $casebam
    #optional
    $star_index


  **Used parameters**

  .. code-block:: bash

    # config.sh
    $outdir
    $ncores
    # asevent_config.sh
    $use_bam_input_files
    $read_length



|tool|
