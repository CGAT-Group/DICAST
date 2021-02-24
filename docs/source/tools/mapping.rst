Mapping tools
=============

Index building takes time
Use mapping config file

=========================  ===============  ==========================
Tool                       Git			   			Comment
=========================  ===============  ==========================
:doc:`mapping/bbmap`
:doc:`mapping/contextmap`                   ``Requires special input``
:doc:`mapping/crac`
:doc:`mapping/dart`
:doc:`mapping/gsnap`
:doc:`mapping/hisat`
:doc:`mapping/mapsplice`
:doc:`mapping/minimap`
:doc:`mapping/segemehl`
:doc:`mapping/star`
=========================  ===============  ==========================

1. Mapping Input Files
^^^^^^^^^^^^^^^^^^^^^^

.. note::
	The paths assume you are using our suggested :doc:`../setup/input`. For example input files see :doc:`examples`.

If not stated otherwise the mapping tools require the following input files:

* Fastq files for paired end mapping:
	.. code-block:: bash

		input/fastq/*yourFastqFile1*_1.fastq
		input/fastq/*yourFastqFile1*_2.fastq
		input/fastq/*yourFastqFile2*_1.fastq
		input/fastq/*yourFastqFile2*_2.fastq
		. . .

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
