General parameters
==================
Found in file: :guilabel:`scripts/config.sh`

Basic parameters
^^^^^^^^^^^^^^^^

ncores
   | Number of cores or threads that each tool will use. Note when using a snakemake pipeline: the resulting number of cores used is a result of multiplication of ncores and snakemake -j parameter. 
   | Default: ``16``

workdir
   | Name of the base directory inside the Docker. 
   | Default: ``/MOUNT``

outdir
   | Name of the output directory; should be named after the specific tool that was used (use the ``$tool`` variable for that).
   | Default: ``$workdir/output/${tool:-unspecific}-output``

read_length
   | Length of the reads inside the fastq files.
   | Default: ``76``

Input Directories
^^^^^^^^^^^^^^^^^

inputdir
   | Base input directory.
   | Default: ``$workdir/input``
   
controlfolder
   | Directory for all needed input files when no differential comparison. Directory for control sample input files when running differential AS event detection.
   | Default: ``$inputdir/controldir``

casefolder
   | Directory only for case sample input files in case of differential AS event detection.
   | Default: ``$inputdir/casedir``

fastqdir
   | Directory for fastq files.
   | Default: ``$controlfolder/fastqdir``

bamdir
   | Directory for bam files.
   | Default: ``$controlfolder/bamdir``

samdir
   | Directory for sam files.
   | Default: ``$controlfolder/bamdir``

fastadir
   | Directory for the fasta file (might vary for certain tools, see mapping or AS-specific config files TODO LINK DOCUMENTATION)
   | Default: ``$inputdir``

gtfdir
   | Directory for the gtf file.
   | Default: ``$inputdir``

gffdir
   | Directory for gff file
   | Default: ``$inputdir``

star_index
   | Folder containing a star index built with the ``$gtf`` and ``$fasta`` files (see below), used by: IRFinder, KisSplice, rMATS
   | Default: ``$workdir/index/star_index``


Input Parameters
^^^^^^^^^^^^^^^^

.. _input_parameters:

fastaname
   | Name of the genome reference file (fasta format) inside ``$fastadir``.
   | Example: ``Homo_sapiens.GRCh38.dna.primary_assembly.fa``

gtfname
   | Name of gtf reference file inside ``$gffdir``.
   | Example: ``splicing_variants.gtf``

gffname
   | Name of gff reference file inside ``$gffdir``.
   | Example: ``splicing_variants.gff3``

.. note::
	There should be no need to edit ``fasta``, ``gtf`` and ``gff`` since they just combine other parameters. 

fasta
   | Full path to the fasta file.
   | Default: ``${fastadir:-unspecific}/$fastaname``

gtf
   | Full path to the gtf file.
   | Default: ``${gtfdir:-unspecific}/$gtfname``

gff
   | Full path to the gff file.
   | Default: ``${gffdir:-unspecific}/$gffname``


ASimulatoR Parameters
^^^^^^^^^^^^^^^^^^^^^

asimulator_inputdir
   | Full path to the input directory used by ASimulatoR.
   | Example: ``/dockers/src/ASimulatoR/in``

asimulator_outputdir
   | Full path to the output directory used by ASimulatoR.
   | Example: ``/dockers/src/ASimulatoR/out``

