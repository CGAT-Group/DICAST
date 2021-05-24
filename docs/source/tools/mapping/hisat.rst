.. Links

.. _manual: http://daehwankimlab.github.io/hisat2/manual/
.. |tool| replace:: HISAT


HISAT 2
=======

.. sidebar:: |tool| Factsheet

  =============  =================
  **Toolname:**  *hisat*
  **Version:**   *2.0.0*
  =============  =================

  **Required files:**

  .. code-block:: bash

    # config.sh
    $fastq
    # only for indexing
    $fasta
    $gtf



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

$fasta
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
	If there is no index it will be automatically built. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.

	.. code-block:: bash

		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh:
		# 	indexdir=hisat_index
		#	  indexname=$fasta_index
		# $fasta to make sure we have the right index for the used fasta file

		index/*your $indexdir variable*/*your $indexname variable*.1.ht2
		index/*your $indexdir variable*/*your $indexname variable*.2.ht2
		index/*your $indexdir variable*/*your $indexname variable*.3.ht2
		. . .

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. The variables can be changed in
:guilabel:`scripts/config.sh` and :guilabel:`scripts/mapping_config.sh`
If you want to specify your analysis with different parameters you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.

  --x
    Basename of the index for the reference genome.

    .. code-block:: bash

      --x $indexname

  -1
    Fastq filename of paired end read 1.

    .. code-block:: bash

      -1 *yourFastqFile1_*1.fastq

  -2
    Fastq filename of paired end read 2.

    .. code-block:: bash

      -2 *yourFastqFile1_*2.fastq

  -S
    The path to the output file in sam format.
    For differential analysis the output will be separated into case and control folder based on the basefolder of the according fastq file.

    .. code-block:: bash

      -S $outdir/controldir/*yourFastqFile1_*.sam

  --known-splicesite-infile
  	Provide a list of known splice sites.

    .. code-block:: bash

      --known-splicesite-infile $indexdir/$indexname/splicesites.txt

  -q
    Activate quiet mode so only error messages are printed.

3. Other comments:
^^^^^^^^^^^^^^^^^^

.. note::

	The index building will also compute a file :guilabel:`splicesites.txt` with known splice sites based on the given gtf file.


4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `Fast and SNP-tolerant detection of complex variants and splicing in short reads <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2844994/>`_
