
.. _snakemake: https://snakemake.readthedocs.io/en/stable/


Install Snakemake
=================

``snakemake`` is the pipe-lining software that enables your DICAST runs. You can set up snakemake in a conda environment.

If you have never worked with conda before you might want to get conda first: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

``Mamba`` is a faster implementation of conda and while optional, is a recommended installation step to support your conda: https://mamba.readthedocs.io/en/latest/installation.html

.. prompt:: bash

	# create conda environment from .yml file with snakemake in it.
	mamba env create -f scripts/snakemake/dicast-snakemake.yml

	# if you want to use DICAST, activate the "dicast-snakemake" environment
	conda activate dicast-snakemake

If you want to learn more about snakemake, you can check out the snakemake documentation: `snakemake`_.


.. toctree::
   :maxdepth: 1
