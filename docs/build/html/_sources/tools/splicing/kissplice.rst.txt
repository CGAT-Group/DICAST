

.. Links

.. _manual: *not available*
.. |tool| replace:: Kissplice

Kissplice
=========

.. note::

  |tool| can be used to calculate differential splicing as well as only alternative-splicing events.
  If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.
  Otherwise set ``differential=0``.

.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *kissplice*
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
    #optional:
    $star_index

  **Used parameters**

  .. code-block:: bash

    # config.sh
    $outdir
    $ncores
    # asevent_config.sh
    $read_length
    $differential

|tool|
