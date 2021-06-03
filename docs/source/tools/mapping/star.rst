.. Links

.. _manual: 
.. |tool| replace:: Star

Star
====



.. sidebar:: |tool| Factsheet

	=============  =========================
	**Toolname:**  *star*
	**Version:**   
	=============  =========================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		$gtf
		# Only for indexing:
		$fasta


Segemehl

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

$gtf
	Gene annotation file in GTF format.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file name

		input/*yourGTFfile*.gtf
    

Optional: Index
	If there is no index it will be automatically built with Segemehl. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh:
		# 	indexdir=|tool|_index
		#	indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

			
	--sjdbGTFfile
		Gene annotation file in GTF format.

		.. code-block:: bash

			-sjdbGTFfile $gtf
		
	--readFilesIn
		Fastq filename of paired end read 1 and pairedend read 2.

		.. code-block:: bash
		
			--readFilesIn *yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq
			
	--genomeDir		
		Basename of Star index files.

   		
		.. code-block:: bash

   			--genomeDir $indexdir/$indexname
		

	--outFileNamePrefix
		Prefix of the output folder.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq files.

		.. code-block:: bash

			--outFileNamePrefix $outdir/controldir/*yourFastqFile1_*|tool|

	--runTreadN
		Set number of threads to be used during the computation

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			# 	$ncores=4

			--runTreadN $ncores
	
	--twopassMode
	
		.. code-block:: bash
			--twopassMode Basic
	
	--outSAMstrandField
		
		.. code-block:: bash
			--outSAMstrandField intronMotif
	
	--outSAMattributes 
		
		.. code-block:: bash
			NH HI AS nM NM XS



3. Other comments:
^^^^^^^^^^^^^^^^^^



4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: 
