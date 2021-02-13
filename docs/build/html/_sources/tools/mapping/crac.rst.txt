.. Links

.. _manual: http://crac.gforge.inria.fr/
.. |tool| replace:: Crac


Crac
====

.. sidebar:: |tool| Factsheet

	=============  =================
	**Toolname:**  *crac*
	**Version:**   *2.4.0*
	=============  =================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		# Only for indexing:
		$fasta



1. Input Files
^^^^^^^^^^^^^^

The following files are required to run |tool|.

.. note::
	The filepaths assume you are using our :doc:`folder structure</setup/input>`.
	Make sure to read the comments (#) as well.

$fastq
	Fastq files for paired end mapping. The path variable can be found in :guilabel:`scripts/config.sh`.
	If you want to perform a non differential analysis put your files inside the :guilabel:`controldir/fastqir/` directory.
	If you want to compare two different experimental settings (case and control), put your fastq files in the :guilabel:`controldir/fastqir/`
	and :guilabel:`casedir/fastqir/` folders respectively.

	.. code-block:: bash

		# Fastq file paths
		# Replace the text between the stars *...* with your file names

		input/controldir/fastqir/*yourFastqFile1_*1.fastq
		input/controldir/fastqir/*yourFastqFile1_*2.fastq
		input/controldir/fastqir/*yourFastqFile2_*1.fastq
		input/controldir/fastqir/*yourFastqFile2_*2.fastq
		. . .

$fasta:
	The name of the reference fasta file to be used (only used for building the index). The path variable can be found in :guilabel:`scripts/config.sh`.

	.. code-block:: bash

		# Fasta files paths
		# Replace the text between the stars *...* with your file name

		input/*yourFastaFile*.fa


Optional: Index
	If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh:
		# 	indexdir=crac_index
		#	indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*.conf
		index/*your $indexdir variable*/*your $indexname variable*.ssa


2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	-i
		Basename of the index-folder/file.

		.. code-block:: bash

			-i $indexdir/$indexname

	-k
		Number of k-mers to be used. 22 is the recommended number for human genome.

		.. code-block:: bash

			-k 22

	-r
		File paths to reads in fastq format. One pair of fastq files for paired-end mapping.

		.. code-block:: bash

			-reads *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq

	-o
		The path to the output directory.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq file.

		.. code-block:: bash

			-o $outdir/controldir/*yourFastqFile1_*crac.sam

	--nb-threads
		Set number of threads to be used during the computation

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			# 	$ncores=4

			--nb-threads $ncores

	--detailed-sam
		Return a detailed sam file as output.

	--stranded
		Reads are from a strand specific RNA-seq protocol



3. Other comments:
^^^^^^^^^^^^^^^^^^



4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `CRAC: an integrated approach to the analysis of RNA-seq reads <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-3-r30>`_
