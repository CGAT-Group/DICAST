Welcome to DICAST!
===================


The DICAST docker collection was initially designed to benchmark alternative splicing tools based on simulated "ground truth" reads produced by ASimulator. Here we provide a pipeline for running several mapping and alternative/differential splicing tools and evaluate and compare the results.  The containers are however not only suitable for simulated data. You can use it for real data too.

.. warning::
	We tried our best to unify the input that is required for all tools. This did not work for all tools. When a tool requires custom input you will see a warning like this on the concerning documentation page.

To provide a fair baseline while maintaining easy usability, per default we run the tools with their default variables. If you feel like this is not doing your tool proper justice please contact us. The default parameters can be changed by editing the ENTRYPOINT.sh scripts of each tool.

The tools included here are the most widely used and well maintained tools in the field. If you would like to include your tool in the pipeline please let us know.
We hope this collection can be a starting point for future benchmarking approaches and quality control.


What is DICAST?
===============

DICAST is a collection of tools for analysing RNA-Seq data. For easy installation and maintenance we provide docker container for every integrated tool.

DICAST can be run as complete pipline, starting with simulating RNA-Seq data with ASimulator, mapping the reads to a fasta reference, get information about alternative splicing with one or multiple tools and finally visualize and compare the results from different tools with DICAST unify. Alternatively you can analyse self-created experimental or simulated RNA-Seq data with the same tools.

When should I use it?
=====================

* If you want to benchmark or compare different mapping and splicing tools
* If you want to analyse your data with one or more tools included in DICAST




.. toctree::
	:maxdepth: 3
	:hidden:
	:caption: Get started

	get-started/setup
	get-started/run
	get-started/advanced


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
