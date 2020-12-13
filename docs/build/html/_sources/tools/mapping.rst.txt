Mapping tools
=============

Index building takes time
Use mapping config file

===================  ===============  ==========================
Tool				 Git			  Comment
===================  ===============  ==========================
:ref:`bbmap`                 
:ref:`contextmap 2`            		  ``Requires special input``
:ref:`crac`                 
:ref:`dart`               
:ref:`gsnap`              
:ref:`hisat 2`            
:ref:`mapsplice 2`        
:ref:`minimap 2`          
:ref:`segemehl`           
:ref:`star`               
===================  ===============  ==========================

1. Mapping Input Files
^^^^^^^^^^^^^^^^^^^^^^

.. note::
	The paths assume you are using our suggested :ref:`folder structure`. For example input files see :ref:`examples`.
	
If not stated otherwise the mapping tools require the following input files:

* Fastq files for paired end mapping:
	.. code-block:: bash
	
		input/fastq/*yourFastqFile1*_1.fastq
		input/fastq/*yourFastqFile1*_2.fastq
		input/fastq/*yourFastqFile2*_1.fastq
		input/fastq/*yourFastqFile2*_2.fastq
		...
* Fasta Reference:
	* input/*yourFastaFile*.fa
* Optional: Index (if there is no index it will be built)
	* index/*yourIndexBaseName*

2. Mapping Config Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^




.. toctree::
   :maxdepth: 2
   :hidden:
   
   mapping/bbmap
   mapping/contextmap
   mapping/crac
   mapping/dart
   mapping/gsnap
   mapping/hisat
   mapping/mapsplice
   mapping/minimap
   mapping/segemehl
   mapping/star
