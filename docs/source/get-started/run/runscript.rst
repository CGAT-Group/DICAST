Use multiple tools at once
==========================

In this section we will explain how to use DICAST to run a whole pipeline.

Make sure you followed the steps described in the setup section carefully.
Before getting started make sure to activate the snakemake conda environment:

.. prompt:: bash

  conda activate snakemake

There are two possibilities to run DICAST as pipeline:

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

  $sequanix -s scripts/snakemake/Snakefile -w ./ -c scripts/snakemake/snakemake_config.yaml


Run pipeline without GUI
^^^^^^^^^^^^^^^^^^^^^^^^^
.. prompt:: bash $,# auto

  See snakemake -h
  E.g.:
  $snakemake -j 2 -d /nfs/proj/AS_dockers/ --configfile snakemake_config.yaml

Troubleshooting
^^^^^^^^^^^^^^^

- Check Logger tab if something isn't working
- Enable tooltips (Option -> Preferences -> Show tooltips) for explanations of the different config parameters
- Snakemake creates a run report. If you want to use the build-in Open Report button, you also have to set the parameter under Option -> Preferences -> HTML page to open as a report. The file name has to be the same .
- Snakemake checks if the output file name is the same as from the run before, if that's the case it won't run again, so be sure to enter a new output file name or delete the old output!


.. toctree::
   :maxdepth: 2
