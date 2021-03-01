

.. Links

.. _manual: *not available*
.. |tool| replace:: Majiq

Majiq
=====

.. note::

  |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
  If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.
  Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *majiq*
  **Version:**   *v*
  **License**    *L*
  =============  =================

  **Required files:**

  .. code-block:: bash

    # config.sh
    $gff
    $controlbam
    # differential analysis only
    $casebam


  **Used parameters**

  .. code-block:: bash

    # config.sh
    $outdir
    $ncores
    # asevent_config.sh
    $differential
    $read_length


|tool|
