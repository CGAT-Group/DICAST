.. Links

.. _manual: https://github.com/hsinnan75/Dart
.. |tool| replace:: Dart
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _license: https://github.com/hsinnan75/Dart/blob/master/LICENSE

Dart
====

.. sidebar:: |tool| Factsheet

	=============  =================
	**Toolname:**  *dart*
	**Version:**   *v1.4.0*
	=============  =================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		$fasta


Dart is a RNA-Seq transcript aligner using a divide and conquer strategy. For more information look up the tools 'manual'_.


1. Input Files
^^^^^^^^^^^^^^
The following files are required to run |tool|.

.. note::
	The filepaths assume you are using our :doc:`folder structure</setup/input>`.
	Make sure to read the comments (#) as well.

$fastq
	Fastq files for paired end mapping. The path variable can be found in :guilabel:`scripts/config.sh`.

	.. code-block:: bash

		# Fastq file paths
		# Replace the text between the stars *...* with your file names

		input/fastq/*yourFastqFile1_*1.fastq
		input/fastq/*yourFastqFile1_*2.fastq
		input/fastq/*yourFastqFile2_*1.fastq
		input/fastq/*yourFastqFile2_*2.fastq
		. . .

$fasta:
	The name of the reference fasta file to be used. The path variable can be found in :guilabel:`scripts/config.sh`.

	.. code-block:: bash

		# Fasta files paths
		# Replace the text between the stars *...* with your file name

		input/*yourFastaFile*.fa

Optional: Index
	|tool| requires a bwt or bwa based index. If there is no index it will be automatically built with our ENTRYPOINT.sh script. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh:
		# 	indexdir=dart_index
		#	indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*.amb
		index/*your $indexdir variable*/*your $indexname variable*.ann
		index/*your $indexdir variable*/*your $indexname variable*.bwt
		index/*your $indexdir variable*/*your $indexname variable*.pac
		index/*your $indexdir variable*/*your $indexname variable*.sa

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	-i
		Prefix of the index to be used.

		.. code-block:: bash

			-i $indexdir/$indexname

	-f
		Fastq filename of paired end read 1.

		.. code-block:: bash

			-f *yourFastqFile1_*1.fastq

	-f2
		Fastq filename of paired end read 2.

		.. code-block:: bash

			-f2 *yourFastqFile1_*2.fastq

	-o
		Output file in sam format.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq file.

		.. code-block:: bash

			-o $outdir/*yourFastqFile1_*dart.sam

3. Other comments:
^^^^^^^^^^^^^^^^^^

Make sure to use a bwt/bwa based index. Other than that, |tool| has no special requirements.

4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `DART: a fast and accurate RNA-seq mapper with a partitioning strategy<https://academic.oup.com/bioinformatics/article/34/2/190/4104410>`_
