Welcome to DICAST!
===================


The DICAST pipeline was initially designed to benchmark alternative splicing (AS) tools based on simulated "ground truth" reads produced by ASimulator. Here we provide a pipeline for running several mapping and AS event detection tools and evaluate and compare the results.  DICAST however is not only suitable for simulated data, but you can also use it for real data.

We hope for this to be an open collaborative effort to represent and benchmark your tools. Please feel free to reach out to us, should your tool already be here and you have some edits to suggest. Also please reach out to us, if your Alternative splicing tool isn't here and you'd like it to be.

To provide a fair baseline while maintaining easy usability, per default we run the tools with their default variables. If you feel like this is not doing your tool justice please contact us. The default parameters can be changed by editing the ENTRYPOINT.sh scripts of each tool.

The tools included here are the most widely used and well maintained among AS event detection tools. If you would like to include your tool in the pipeline please let us know.
We hope this collection can be a starting point for future benchmarking approaches and quality control.

.. figure:: /img/Screenshot_dicast.png


What is DICAST?
===============
.. _Snakemake: https://snakemake.github.io/
.. _Docker: https://www.docker.com/
.. _Dockerhub: https://hub.docker.com/repository/docker/dicastproj/dicast
.. _bioRxiv: https://doi.org/10.1101/2022.01.05.475067

DICAST is a collection of alternative splicing event detection tools for analyzing RNA-Seq data. DICAST runs on `Snakemake`_ pipelines and relies on `Docker`_ based containerization. For easy installation and maintenance, we provide docker containers for every integrated tool at `Dockerhub`_

DICAST can be run as a complete pipeline, starting with simulating RNA-Seq data with ASimulator, mapping the reads to a fasta reference, get information about AS with one or multiple tools and finally visualize and compare the results from different tools with DICAST unify.

Alternatively, you can run one of the same tools as a single docker container without snakemake.

When should I use DICAST?
=========================

* If you want to benchmark or compare different mapping and splicing tools with a genome and an annotation.
* If you want to analyze your data with one or more tools included in DICAST and find out which tools give you the events you're interested in.

How do I cite DICAST?
=====================
DICAST is currently hosted on `bioRxiv`_, at this link `<https://doi.org/10.1101/2022.01.05.475067>`_.



.. toctree::
	:maxdepth: 3
	:hidden:
	:caption: Get started

	get-started/setup
	get-started/config
	get-started/run


.. toctree::
	:maxdepth: 3
	:hidden:
	:caption: Tools

	tools/tools
	tools/mapping
	tools/splicing


.. toctree::
	:maxdepth: 3
	:hidden:
	:caption: Further information

	further-information/example-files
	further-information/faq
	further-information/uninstall-dicast	
	further-information/about
