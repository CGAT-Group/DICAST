
.. _snakemake: https://snakemake.readthedocs.io/en/stable/


Install Snakemake
=================

``snakemake`` is the pipe-lining software that enables your DICAST runs. You can set up snakemake in a conda environment.

If you have never worked with conda before you might want to get conda first: https://conda.io/projects/conda/en/latest/user-guide/install/index.html


.. prompt:: bash

	# create conda environment from .yml file with snakemake in it.
	conda env create -f scripts/snakemake/dicast-snakemake.yml

	# if you want to use DICAST, activate the "dicast-snakemake" environment
	conda activate dicast-snakemake

If you want to learn more about snakemake, you can check out the snakemake documentation: `snakemake`_.


.. toctree::
   :maxdepth: 1
