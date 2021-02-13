.. Links

.. _manual: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/
.. |tool| replace:: BBMap


BBMap
=====

.. sidebar:: |tool| Factsheet

	=============  =================
	**Toolname:**  *bbmap*
	**Version:**   *38.43*
	=============  =================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		$fasta

BBMap uses a multi kmer seed and extend strategy for read mapping.

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
		# 	indexdir=bbmap_index
		#	indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*/ref/genome/
		index/*your $indexdir variable*/*your $indexname variable*/ref/index/

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	-in
		Fastq filename of paird end read 1.

		.. code-block:: bash

			-in *yourFastqFile1_*1.fastq

	-in2
		Fastq filename of paird end read 2.

		.. code-block:: bash

			-in2 *yourFastqFile1_*2.fastq

	-ref
		Reference genome in fasta format.

		.. code-block:: bash

			-ref $fasta

	-path
		Basename of the index-folder/file.

		.. code-block:: bash

			-path $indexdir/$indexname

	-outm
		The path to the **mapped** output file in sam format.
		For differential analysis the output will be separated into case and sontrol folder based on the basefolder of the according fastq file.

		.. code-block:: bash

			-o $outdir/*yourFastqFile1_*bbmap.sam

	-outu
		The path to the **unmapped** output file in sam format.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq file.

		.. code-block:: bash

			-o $outdir/controldir/*yourFastqFile1_*bbmap_unmapped.sam



3. Other comments:
^^^^^^^^^^^^^^^^^^



4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `Long Read RNA-seq Mapper <http://bib.irb.hr/datoteka/773708.Josip_Maric_diplomski.pdf>`_
