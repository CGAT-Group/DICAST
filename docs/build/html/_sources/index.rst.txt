Welcome to DICAST!
===================

The DICAST docker collection was initially designed to benchmark alternative splicing tools based on simulated "ground truth" reads produced by ASimulator. Here we provide a pipeline for running several mapping and alternative/differential splicing tools and evaluate and compare the results.  The containers are however not only suitable for simulated data. You can use it for real data too. 

.. warning::
	We tried our best to unify the input that is required for all tools. This did not work for all tools. When a tool requires custom input you will see a warning like this on the concerning documentation page.

To provide a fair baseline while maintaining easy usability, per default we run the tools with their default variables. If you feel like this is not doing your tool proper justice please contact us. The default parameters can be changed by editing the ENTRYPOINT.sh scripts of each tool. 
To learn more about this, please look at the "Change ENTRYPOINTS" section.

The tools included here are the most widely used and well maintained tools in the field. If you would like to include your tool in the pipeline please let us know.
We hope this collection can be a starting point for future benchmarking approaches and quality control.


.. toctree::
   :maxdepth: 2
   :caption: Table of contents:
   
   setup
   run
   tools
   examples
   faq

Contributions
=============

Funding
=======

Contact
=======

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
