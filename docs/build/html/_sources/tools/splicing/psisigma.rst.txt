
.. Links

.. _manual: *not available*
.. |tool| replace:: psisigma

psisigma
========

.. warning::

	|tool| needs properly named bam folders to work. Please make sure to store the files in the folders corresponding to the following variables in the :guilabel:`/scripts/asevent_config.sh`:
	``renamed_casebam``,``renamed_controlbam`` and ``renamed_star_alignment_files``.

.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *psisigma*
  **Version:**   *v*
  **License**    *L*
  =============  =================

  **Required files:**

  .. code-block:: bash

    # config.sh
    $gtf
    # asevent_config.sh
    $renamed_casebam
    $renamed_controlbam
    $renamed_star_alignment_files


  **Used parameters**

  .. code-block:: bash

    # config.sh
    $outdir


|tool|
