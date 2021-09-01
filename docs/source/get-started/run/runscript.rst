Use multiple tools at once
==========================

In this section we will explain how to use DICAST to run a whole pipeline.

Make sure you followed the steps described in the :doc:`setup <../setup>` section carefully.
Before getting started make sure to activate the snakemake conda environment:

.. prompt:: bash

  conda activate snakemake

.. note::

	Snakemake is set up to run all tools in parallel, meaning if the pipeline is run unrestricted, it will use all available cores. Be sure to set the number of cores to limit the resources available to Snakemake (explained below).



There are two possibilities to run the DICAST pipeline:

  - With graphical user interface (recommended)
  - With command line

Run pipeline with GUI
^^^^^^^^^^^^^^^^^^^^^^

.. prompt:: bash $,# auto

  sequanix [OPTIONS]
    -w WKDIR, --working-directory WKDIR
                          Set working directory
    -s SNAKEFILE, --snakefile SNAKEFILE
                          A valid Snakefile
    -c CONFIGFILE, --configfile CONFIGFILE

  E.g.:

  $ sequanix -s scripts/snakemake/Snakefile -w ./ -c scripts/snakemake/snakemake_config.yaml

This command will open the GUI from which the pipeline can be started. Make sure the *Generic pipelines* tab is selected. 
Also be sure to set the desired number of cores to be used by Snakemake by going to *Option* -> *Snakemake Options* -> *Local* or press Ctrl + O.

For more information, see the `sequanix documentation <https://sequana.readthedocs.io/en/master/sequanix.html>`_. 

TODO: More detailed explanation once we have the final GUI (@Zaka)

Run pipeline without GUI
^^^^^^^^^^^^^^^^^^^^^^^^^
.. prompt:: bash $,# auto

  See snakemake -h
  E.g.:
  $ snakemake -j 2 -d /nfs/proj/AS_dockers/ --configfile snakemake_config.yaml

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

Troubleshooting
^^^^^^^^^^^^^^^

- Check Logger tab if something isn't working
- If a run was canceled/exited unexpectetly and the directory is still locked, try running the same snakemake command again with ``--unlock`` (or use the Unlock button in the GUI) or remove the files from ``working_directory/.snakemake/locks/``.
- Snakemake itself creates some logging files, they can be found in ``working_directory/.snakemake/log/``
- Enable tooltips (Option -> Preferences -> Show tooltips) for explanations of the different config parameters
- Snakemake creates a run report. If you want to use the build-in Open Report button, you also have to set the parameter under Option -> Preferences -> HTML page to open as a report. The file name has to be the same.

.. toctree::
   :maxdepth: 2
