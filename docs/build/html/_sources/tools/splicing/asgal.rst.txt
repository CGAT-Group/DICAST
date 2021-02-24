.. Links

.. _manual: *not available*
.. |tool| replace:: Asgal

Asgal
=====

.. warning::

	Asgal requires the variables ``$fastqpair1suffix`` and ``$fastqpair2suffix`` to be set in the :guilabel:`scripts/asevent_config.sh` file.

.. sidebar:: |tool| Factsheet

	=============  =================
	**Toolname:**  *contextmap*
	**Version:**   *v1.1.1*
	=============  =================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastqdir/*$fastqpair1suffix
		$fastqdir/*$fastqpair2suffix
		$fasta
		$transcript
		$gtf

|tool|

1. Input Files
^^^^^^^^^^^^^^

The following files are required to run |tool|.

.. note::
	The filepaths assume you are using our :doc:`folder structure</setup/input>`.
	Make sure to read the comments (#) as well.

$fastqdir/$fastqpair1suffix
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


$fastqdir/$fastqpair2suffix
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
	Gene annotation file in GTF format. The path variable can be found in :guilabel:`scripts/config.sh`.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file names

		input/*yourGTFfile*.gtf

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/asevent_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	--multi


		.. code-block:: bash

			--multi

	-g
		Fasta reference file.

		.. code-block:: bash

			-g $fasta

	-a
		Gtf annotation file.

		.. code-block:: bash

			-a $gtf

	-t
		Gene transcripts in fasta format

		.. code-block:: bash

			-t $transcript

	-s
		Fastq filename of paired end read 1.

		.. code-block:: bash

			-f *yourFastqFile1_*1.fastq

	-s2
		Fastq filename of paired end read 2.

		.. code-block:: bash

			-f2 *yourFastqFile1_*2.fastq

	-o
		The path to the output directory for the accourding fastq file pair. The file will be named after the fastq file basename.

		.. code-block:: bash

			-o $outdir/*yourFastqFile1*-output

	-@
		Set number of threads to be used during the computation.

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			#    $ncores=4

			-@ $ncores

3. Other comments:
^^^^^^^^^^^^^^^^^^


4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: 
