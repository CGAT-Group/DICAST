.. Links

.. _manual: http://www.netlab.uky.edu/p/bioinfo/MapSplice2UserGuide
.. |tool| replace:: Mapsplice
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _latest conda version: https://anaconda.org/bioconda/mapsplice

Mapsplice 2
===========

.. warning::

	Mapsplice requires a :ref:`special input<mapsplice_input>`: ``$mapsplice_fastadir_index`` and ``mapsplice_fastadir_mapping``.


.. sidebar:: |tool| Factsheet

	=============  =========================
	**Toolname:**  *mapsplice*
	**Version:**   `latest conda version`_
	=============  =========================

	**Required files:**

	.. code-block:: bash

		# config.sh
		$fastq
		$gtf
		# mapping_config.sh
		$mapsplice_fastadir_index
		$mapsplice_fastadir_mapping



Mapsplice

.. _mapsplice_input:

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

$gtf
	Gene annotation file in GTF format.

	.. code-block:: bash

		# Replace the text between the stars *...* with your file names

		input/*yourGTFfile*.gtf

$mapsplice_fastadir_index:
	|tool| requires chromosome-wise fasta files as reference for the **indexing**. The path variable can be found in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Chromosome-wise fasta files paths
		# Replace the text between the stars *...* with your file names

		input/fasta_chromosomes/*yourChromosome1*.fa
		input/fasta_chromosomes/*yourChromosome2*.fa
		input/fasta_chromosomes/*yourChromosome3*.fa
		. . .

$mapsplice_fastadir_mapping
	|tool| requires a special formatted fasta file for the **mapping** step.
	The path variable can be found in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Special formatted mapsplice fasta file.
		# Replace the text between the stars *...* with your file names

		input/fasta_mapsplice/*yourFastaForMapsplice*.fa

Optional: Index
	If there is no index it will be automatically built with bowtie1. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh:
		# 	indexdir = mapsplice_index
		#	  indexname = $fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*.1.ebwt
		index/*your $indexdir variable*/*your $indexname variable*.2.ebwt
		index/*your $indexdir variable*/*your $indexname variable*.3.ebwt
		index/*your $indexdir variable*/*your $indexname variable*.4.ebwt
		index/*your $indexdir variable*/*your $indexname variable*.rev.1.ebwt
		index/*your $indexdir variable*/*your $indexname variable*.rev.2.ebwt


2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

	-c
		Mapsplice reference genome in fasta format.

		.. code-block:: bash

			-c $inputdir/$mapsplice_fastadir_mapping

	-x
		Basename of the bowtie1 index.

		.. code-block:: bash

			-x $indexdir/$indexname

	--gene-gtf
		The path to the gene annotation file in GTF format for annotation of fusion junctions.

		.. code-block:: bash

			--gene-gtf $gtf

	-o
		The path to the output file in sam format.
		For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq files.

		.. code-block:: bash

			-o $outdir/controldir/*yourFastqFile1_*mapsplice.sam

	-p
		Set number of threads to be used during the computation.

		.. code-block:: bash

			# If you use our default parameters and folder structure:
			# 	$ncores=4

			-p $ncores

	-1
		Fastq filename of paired end read 1.

		.. code-block:: bash

			-1 *yourFastqFile1_*1.fastq

	-2
		Fastq filename of paired end read 2.

		.. code-block:: bash

			-2 *yourFastqFile1_*2.fastq

3. Other comments:
^^^^^^^^^^^^^^^^^^

TODO: Add info on .fa file

4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `MapSplice: Accurate mapping of RNA-seq reads for splice junction discovery <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2952873/>`_
