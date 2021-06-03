Mapping tools
=============

.. note::

	Most mapping tools need an index file of the reference genome for mapping. The computation of these files can take a long time.
	Our mapping tool scripts check if there already is an index for the respective tool and only build it, if it is not found.
	If you face any index related errors, either set the parameter ``$recompute_index=True`` or delete the old index to recalculate it.

* Index building takes time
* Use mapping config file

=========================  ===============  ==========================
Tool                       Git			   			Comment
=========================  ===============  ==========================
:doc:`mapping/bbmap`
:doc:`mapping/contextmap`                   ``Requires special input``
:doc:`mapping/crac`
:doc:`mapping/dart`
:doc:`mapping/gsnap`
:doc:`mapping/hisat`
:doc:`mapping/mapsplice`                    ``Requires special input``
:doc:`mapping/minimap`                      ``Requires special input``
:doc:`mapping/segemehl`
:doc:`mapping/star`
=========================  ===============  ==========================

1. Mapping Input Files
^^^^^^^^^^^^^^^^^^^^^^

.. tip::
	The paths assume you are using our suggested :doc:`input structure <../get-started/setup/input>`.
	Example input files you can find in our :doc:`examples section<../further-information/example-files>`.

If not stated otherwise the mapping tools require the following input files:

* Fastq files for paired end mapping:
	.. code-block:: bash

		input/fastq/*yourFastqFile1*_1.fastq
		input/fastq/*yourFastqFile1*_2.fastq
		input/fastq/*yourFastqFile2*_1.fastq
		input/fastq/*yourFastqFile2*_2.fastq
		. . .

* Fasta Reference:
	.. code-block:: bash

		input/*yourFastaFile*.fa

* Gtf File:
	.. code-block:: bash

		input/*yourGtf*.gtf

* Optional: Index
  Tool specific index file(s).
	If you want to provide your own index please make sure it is in the correct format and file names.
	If no index file is found in the index folder it will be built the first time you run the tool.
	.. code-block:: bash

		index/*yourIndexBaseName*

2. Special input files
^^^^^^^^^^^^^^^^^^^^^^

Some mapping tools require special input files as indicated in the table above.
For more information please consult the respective documentation page.

3. Parameters
^^^^^^^^^^^^^^




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
