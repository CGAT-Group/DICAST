Splicing tools
==============

Differential splicing vs alternative splicing
Use asevent config file

Input Files
^^^^^^^^^^^^^^^^^^^^^^

.. note::

	Not all tools need the same input files.
	Check out the respective tool documentation sites to find out which files are needed to run the tool.

$fastqdir/.\*$fastqpair1suffix
	Fastq files for pair 1 fastq files stored in ``$fastqdir``, identified by the suffix ``$fastqpair1suffix``.
	The path variables can be found in :guilabel:`scripts/config.sh` and :guilabel:`scripts/asevent_config.sh`.

	.. code-block:: bash

		# Fastq file paths
		# Assumed variable settings:
		#    $fastqdir=input/fastq   ## in config.sh
		#    $fastqpair1suffix="_1.fastq"   ## in asevent_config.sh
		# Replace the text between the stars *...* with your file names

		input/fastq/*yourFastqFile1*_1.fastq
		input/fastq/*yourFastqFile2*_1.fastq
		. . .


$fastqdir/.\*$fastqpair2suffix
	Fastq files for pair 1 fastq files stored in ``$fastqdir``, identified by the suffix ``$fastqpair2suffix``.
	The path variables can be found in :guilabel:`scripts/config.sh` and :guilabel:`scripts/asevent_config.sh`.

	.. code-block:: bash

		# Fastq file paths
		# Assumed variable settings:
		#    $fastqdir=input/fastq   ## in config.sh
		#    $fastqpair2suffix="_2.fastq"   ## in asevent_config.sh
		# Replace the text between the stars *...* with your file names

		input/fastq/*yourFastqFile1*_2.fastq
		input/fastq/*yourFastqFile2*_2.fastq
		. . .

$controlbam

casebam

$fasta:
	The name of the reference fasta file. The path variable can be found in :guilabel:`scripts/config.sh`.

	.. code-block:: bash

		# Fasta files paths
		# Replace the text between the stars *...* with your file name

		input/*yourFastaFile*.fa

$transcript
	The name of the fasta file for gene transcripts. The path variable can be found in :guilabel:`scripts/asevent_config.sh`.

	.. code-block:: bash

		# Assumed variable settings:
		#    $inputdir=input   ## in config.sh

		index/*yourTranscriptFasta*.fasta

$gtf
	Gene annotation file in GTF format.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file name

		input/*yourGTFfile*.gtf

$gff
	Gene annotation file in GFF format.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file name

		input/*yourGFFfile*.gff

Input Parameters
^^^^^^^^^^^^^^^^^^^^^^

.. note::

	Not all tools use all the parameters
	Check out the respective tool documentation sites to find out which parameters are used by a tool.


$ncores
	Number of threads to be used during the computation. Set to ``ncores=4`` per default.

$differential
  Some tools can be used to calculate differential splicing as well as only alternative-splicing events.
  If you want to perform differential analysis set ``differential=1`` in the :guilabel:`/scripts/asevent_config.sh` config file.
  It is set to ``differential=0`` per default, which indicates AS event detection.

$read_length
  Read length of reads in your fastq files.





.. toctree::
   :maxdepth: 1
   :hidden:

   splicing/asgal
   splicing/aspli
   splicing/cash
   splicing/dartsbht
   splicing/dexseq
   splicing/dsplicetype
   splicing/eventpointer
   splicing/irfinder
   splicing/jum
   splicing/juncbase
   splicing/kissplice
   splicing/leafcutter
   splicing/majiq
   splicing/miso
   splicing/psisigma
   splicing/rmats
   splicing/sajr
   splicing/sgseq
   splicing/spladder
   splicing/whippet
