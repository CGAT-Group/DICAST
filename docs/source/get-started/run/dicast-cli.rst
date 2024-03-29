Run DICAST via Command Line Interface 
======================================

In this section we will explain how to use DICAST to run a whole pipeline on the terminal alone.

Make sure you followed the steps described in the :doc:`setup <../setup>` section carefully.

Change config.sh according to your run (see :doc:`How to change your config.sh file <../config>`)

Before getting started make sure to activate the snakemake conda environment:

.. prompt:: bash

  conda activate dicast-snakemake

.. note::

	Snakemake is set up to run all tools in parallel, meaning if the pipeline is run unrestricted, it will use all available cores. Be sure to set the number of cores to limit the resources available to Snakemake (explained below).



To run the snakemake pipeline:
 - Go to /path/to/DICAST/scripts/snakemake
 - Edit snakemake_config.yaml to list the tools you want to run in the corresponding lines.
 - Use a snakemake command. For example:

.. prompt:: bash $,# auto

  See snakemake -h
  E.g.:
  $ snakemake -j 2 -d /opt/DICAST/ -s /path/to/DICAST/scripts/snakemake/Snakefile-cli --configfile /path/to/DICAST/scripts/snakemake/snakemake-cli_config.yaml

Important arguments:

=====================================  ==========================================================================================
Argument                               Explanation
=====================================  ==========================================================================================
``-j``, ``--cores <n>``                Set the number of cores used by the pipeline.
``-s``, ``--snakefile <file>``         Set the path to the Snakemake file.
``-d``, ``--directory <path>``         Set the path to the working directory (containing the input, scripts folder etc.)
``--configfile <file>``                Set the path to the configuration file which specifies what part of the pipeline to run.
=====================================  ==========================================================================================

For more information, see the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_.

.. warning::

  This feature has been depreciated and needs another Snakefile, without the ``Possible_overwrite_acknowledge`` rule. We keep this documentation here, for future support and to show you how the tool works under the hood.

Troubleshooting
^^^^^^^^^^^^^^^

- Check log files under output/<tool>-output/logs/
- If a run was canceled/exited unexpectetly and the directory is still locked, try running the same snakemake command again with ``--unlock`` or remove the files from ``working_directory/.snakemake/locks/``
- Snakemake itself creates some logging files, they can be found in ``working_directory/.snakemake/log/``

.. toctree::
   :maxdepth: 2
