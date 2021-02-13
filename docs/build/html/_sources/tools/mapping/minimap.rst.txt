.. Links

.. _manual: https://lh3.github.io/minimap2/minimap2.html
.. |tool| replace:: Minimap
.. _latest conda version: https://anaconda.org/bioconda/minimap

Minimap 2
=========



.. sidebar:: |tool| Factsheet

	=============  =========================
	**Toolname:**  *minimap*
	**Version:**   `latest conda version`_
	=============  =========================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		# Only for indexing:
		$fasta


Minimap

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
	If there is no index it will be automatically built with bowtie1. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh:
		# 	indexdir=minimap_index
		#		indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	-a
		Generate CIGAR and provide output in sam format (PAF is default for |tool|).

	-o
		The path to the output file in sam format.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq files.

		.. code-block:: bash

			-o $outdir/controldir/*yourFastqFile1_*minimap.sam

	-t
		Set number of threads to be used during the computation

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			# 	$t=4

			-t $ncores

	reads
		After all other options call reads in fastq format. One pair of fastq files for paired-end mapping.

		.. code-block:: bash

			*yourFastqFile1_*1.fastq *yourFastqFile1_*2.fastq

3. Other comments:
^^^^^^^^^^^^^^^^^^

TODO: Add info on .fa file

4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `Minimap2: pairwise alignment for nucleotide sequences <https://doi.org/10.1093/bioinformatics/bty191>`_
