Splicing tools
==============

.. warning::
	Currently **only alternative splicing event detection is supported**. Differential splicing tools are coming soon. The differential splicing function of tools which are able to compute both alternative and differential splicing the differential mode is still in beta.

For splicing tools we differentiate between alternative and differential splicing tools. Some tools are able to compute both.
Differential splicing tools compute alternative splicing for two conditions (e.g. case and control) and the files should be separated as indicated by our input folder structure. For alternative splicing analysis "control" is the default.

Splicing Input Files
^^^^^^^^^^^^^^^^^^^^^^^

.. tip::
  The paths assume you are using our suggested :doc:`input structure <../get-started/setup/input>`.
  Example input files you can find in our :doc:`examples section<../further-information/example-files>`.

You can find the required input files in the tool-specific documentation.

.. _fastqSplicing:

fastq
	Fastq files for pair 1 and 2 fastq files stored in ``$fastqdir``, identified by the suffix ``$fastqpair1suffix`` and ``$fastqpair2suffix`` respectively. Not all splicing tools work with fastq files.
	The path variables can be found in :guilabel:`scripts/config.sh` and :guilabel:`scripts/asevent_config.sh`. For differential splicing the files need to be separated in ``controldir`` and ``casedir``

	.. code-block:: bash

		# Fastq file paths
		# Assumed variable settings:
		#    $fastqdir=input/fastq   ## in config.sh
		#    $fastqpair1suffix="_1.fastq"   ## in asevent_config.sh
		#    $fastqpair2suffix="_2.fastq"   ## in asevent_config.sh
		# Replace the text between the stars *...* with your file names

		input/controldir/fastq/*yourFastqFile1*_1.fastq
		input/controldir/fastq/*yourFastqFile1*_2.fastq
		input/controldir/fastq/*yourFastqFile2*_1.fastq
		input/controldir/fastq/*yourFastqFile2*_2.fastq
		. . .

.. _bamSplicing:

bam
	Bam files created by a mapping tool of your choice. When DICAST is run as a pipeline, these will be created by the selected mapping tool(s).

	.. code-block:: bash

		input/controldir/fastq/*yourFastqFile1*_1.fastq

.. _fastaSplicing:

fasta:
	The name of the reference fasta file. The path variable can be found in :guilabel:`scripts/config.sh`.

	.. code-block:: bash

		# Fasta files paths
		# Replace the text between the stars *...* with your file name

		input/*yourFastaFile*.fa

.. _transcriptSplicing:

transcript
	The name of the fasta file for gene transcripts. The path variable can be found in :guilabel:`scripts/asevent_config.sh`.

	.. code-block:: bash

		# Assumed variable settings:
		#    $inputdir=input   ## in config.sh

		input/*yourTranscriptFasta*.fasta

.. _gtfSplicing:

gtf
	Gene annotation file in GTF format.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file name

		input/*yourGTFfile*.gtf

.. _gffSplicing:

gff
	Gene annotation file in GFF format.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file name

		input/*yourGFFfile*.gff

Parameters
^^^^^^^^^^^^^

To provide a fair baseline while maintaining easy usability, per default we run the tools with their default variables. The default parameters can be changed by editing the ENTRYPOINT.sh scripts of each tool. The variables used by mapping ENTRYPOINT.sh scripts can be set in the ``config.sh`` and ``asevent_config.sh`` files in the ``scripts`` folder. For a usual analysis you should not need to change these parameters.



.. toctree::
   :maxdepth: 1
   :hidden:

   splicing/asgal
   splicing/aspli
   splicing/eventpointer
   splicing/irfinder
   splicing/kissplice
   splicing/majiq
   splicing/sgseq
   splicing/spladder
   splicing/whippet
