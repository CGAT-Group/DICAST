

.. Links

.. _manual: *not available*
.. |tool| replace:: Aspli

Aspli
=====

.. note::

  |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
  If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.
  Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *aspli*
  **Version:**   *v*
  **License**    *L*
  =============  =================

  **Required files:**

  .. code-block:: bash

    # config.sh
    $gtf
    $controlbam
    # differential analysis only
    $casebam

  **Used parameters**

  .. code-block:: bash

    # config.sh
    $differential
    $ncores
    $read_length
    $outdir



|tool|
