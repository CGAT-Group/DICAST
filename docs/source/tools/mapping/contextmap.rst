.. Links

.. _manual: https://www.bio.ifi.lmu.de//files/Software/ContextMap/manual/ContextMap-manual.html
.. |tool| replace:: ContextMap
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

ContextMap 2
============

.. warning::
	
	ContextMap requires a :ref:`special input<1. Input Files>`: ``$contextmap_fastadir``. 
	ALSO: make sure to put the file ``jre-8u241-linux-i586.tar.gz`` in the same folder as the ContextMap Dockerfile (:guilabel:`src/contextmap2/`). You can get the file here: https://www.oracle.com/java/technologies/javase/javase8u211-later-archive-downloads.html

.. sidebar:: |tool| Factsheet
	
	=============  =================
	**Toolname:**  *contextmap*
	**Version:**   *v2_7_9*
	=============  =================
	
	**Required files:**
	
	.. code-block:: bash
			
		# config.sh
		$fastq
		# mapping_config.sh
		$contextmap_fastadir
	
	

ContextMap is a mapping tool for RNA-seq data which uses context information of near-by mapped reads for the read mapping process.


1. Input files
^^^^^^^^^^^^^^
The following files are required to run |tool|. 

.. note::
	The filepaths assume you are using our :ref:`folder structure<Folder Structure>`.
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

$contextmap_fastadir:
	|tool| requires chromosome-wise fasta files as reference. The path variable can be found in :guilabel:`scripts/mapping_config.sh`.
	
	.. code-block:: bash
		
		# Chromosome-wise fasta files paths 
		# Replace the text between the stars *...* with your file names
		
		input/fasta_chromosomes/*yourChromosome1*.fa
		input/fasta_chromosomes/*yourChromosome2*.fa
		input/fasta_chromosomes/*yourChromosome3*.fa
		. . .

Optional: Index 
	If there is no index it will be automatically built with bowtie2. If you want to rebuild the index anyway set ``$recompute_index=true`` in :guilabel:`scripts/mapping_config.sh`.
	Once the 
	
	.. code-block:: bash
		
		# Index files paths
		# Replace the text between the stars *...* with your file names
		# Default variable settings in mapping_config.sh: 
		# 	indexdir=contextmap_index
		#	indexname=$fasta 
		# $fasta to make sure we have the right index for the used fasta file
		
		index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.1.bt2
		index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.2.bt2
		index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.3.bt2
		index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.4.bt2
		index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.rev.1.bt2
		index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.rev.1.bt2
		. . .

2. Default parameters:
^^^^^^^^^^^^^^^^^^^^^^
The following parameters are set in the ENTRYPOINT.sh script in our docker to run |tool|. 
If you want to specify your analysis you will have to change the ENTRYPOINT script.
For further information please consult the |tool| `manual`_.
	
	-reads
		Comma separated list of file paths to reads in fastq format. One pair of fastq files for paired-end mapping
		
		.. code-block:: bash
		
			-reads *yourFastqFile1_*1.fastq,*yourFastqFile1_*2.fastq
			
	-aligner_name
		Used aligner (index tool). We use `bowtie2`_.
		
		.. code-block:: bash
		
			bowtie2
			
	-aligner_bin
		Path to the used aligner. If you use our docker you will not have to wolly about it.
		
		.. code-block:: bash
			
			/home/biodocker/bin/bowtie2
			
	-indexer_bin 	
		Path to the indexing tool of the aligner.
		
		.. code-block:: bash
			
			/home/biodocker/bin/bowtie2-build
		
	-indices
		Comma separated list to your index files base names.
		
		.. code-block:: bash
			
			# e.g. for chromosome 1:
			# If your index files are named as follows:
			# 	index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.1.bt2
			# 	index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.2.bt2
			#	index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.3.bt2
			#	index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*.4.bt2
			# *IndexChromosome1* is:
			#	index/*your $indexdir variable*/*your $indexname variable*/*yourChromosome1*
			
			*IndexChromosome1*,*IndexChromosomes2*,*IndexChromosome3*, . . .
			
			
	-genome
		Directory path with chromosome-wise fasta files.
		
		.. code-block:: bash
		
			# If you use our default parameters and folder structure:
			# 	$contestmap_fastadir=input/fasta_chromosomes
			
			$contextmap_fastadir
			
	-o
		The path to the output directory.
		
		.. code-block:: bash
		
			$outdir/*yourFastqFile1_*contextmap

3. Other comments:
^^^^^^^^^^^^^^^^^^

.. warning::
	The jre file ``jre-8u241-linux-i586.tar.gz`` can't be provided by us due to ORACLE lizense. You will have to download it yourself: https://www.oracle.com/java/technologies/javase/javase8u211-later-archive-downloads.html
	Store the file inside the :guilabel:`src/contextmap2/` folder.

4. Important links:
^^^^^^^^^^^^^^^^^^^
	- |tool| `manual`_
	- |tool| publication: `A context-based approach to identify the most likely mapping for RNA-seq experiments <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-S6-S9>`_
	- `bowtie2`_ (used indexing tool)