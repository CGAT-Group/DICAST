.. _conda: https://docs.conda.io/projects/conda/en/latest/index.html
.. _snakemake: https://snakemake.readthedocs.io/en/stable/


Install Snakemake
=================

.. note::

	If you do not plan on running DICAST as a pipeline but only want to run one tool at a time you do not need snakemake.

You can set up snakemake in a conda environment.
If you have never worked with anaconda before you might want to check out their documentation first: `conda`_.

.. prompt:: bash

	# go to the git Directory
	cd dockers

	# create conda environment from .yml file
	conda env create -f scripts/snakemake/snakemake_conda_env.yml
	
	# if you want to use it: activate environment
	conda activate snakemake

If you want to learn more about snakemake, you can check out the snakemake documentation: `snakemake`_.


.. toctree::
   :maxdepth: 1
