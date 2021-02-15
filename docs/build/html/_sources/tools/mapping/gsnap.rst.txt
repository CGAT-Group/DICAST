.. Links

.. _manual: http://research-pub.gene.com/gmap/
.. |tool| replace:: GSNAP


GSNAP
=====

.. sidebar:: |tool| Factsheet

	=============  =================
	**Toolname:**  *gsnap*
	**Version:**   *38.43*
	=============  =================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		$fasta



1. Input Files
^^^^^^^^^^^^^^

The following files are required to run |tool|.

.. note::
	The filepaths assume you are using our :doc:`folder structure</setup/input>`.
	Make sure to read the comments (#) as well.

$fastq
	Fastq files for paired end mapping. The path variable can be found in :guilabel:`scripts/config.sh`.
	If you want to perform a non differential analysis put your files inside the :guilabel:`controldir/fastqir/` directory
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
		# 	indexdir=gsnap_index
		#	indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*.maps
		index/*your $indexdir variable*/*your $indexname variable*.chromosome
		index/*your $indexdir variable*/*your $indexname variable*.chromosome.iit
		. . .

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	--db
		Genome database (base name of index).

		.. code-block:: bash

			--db $indexname

	--dir
		Genome directory (directory of index).

		.. code-block:: bash

			--dir $indexdir

	--output-file
		The path to the output file in sam format.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq file.

		.. code-block:: bash

			-o $outdir/controldir/*yourFastqFile1_*gsnap.sam

	--format
		Define output format (one of sam, m8)

		.. code-block:: bash

			--format sam

	--nthreads
		Set number of threads to be used during the computation

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			# 	$ncores=4

			--nthreads $ncores

	reads
		After all other options call reads in fastq format. One pair of fastq files for paired-end mapping.

		.. code-block:: bash

			*yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq

3. Other comments:
^^^^^^^^^^^^^^^^^^



4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `Fast and SNP-tolerant detection of complex variants and splicing in short reads <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2844994/>`_
