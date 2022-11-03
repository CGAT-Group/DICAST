Tool core parameters
=======================
Found in file: :guilabel:`scripts/config.sh`

.. note::

  If a parameter is recommended as a default. It's for the snakemake workflow to work smooth. Parameters with this value will be marked with: ``recommended to leave at default``

 .. warning::

   If a parameter exists in the config files but isn't listed in this reference, please don't change the default on this paramenter.


Basic parameters
^^^^^^^^^^^^^^^^

ncores
   | Number of cores or threads that each tool will use. Note when using a snakemake pipeline: the resulting number of cores used is a result of multiplication of ncores and snakemake -j parameter.
   | Default: ``16``

workdir ``recommended to leave at default``
   | Name of the base directory inside the Docker.
   | Default: ``/MOUNT``

outdir ``recommended to leave at default``
   | Name of the output directory; should be named after the specific tool that was used (use the ``$tool`` variable for that).
   | Default: ``$workdir/output/${tool:-unspecific}-output``

read_length
   | Length of the reads inside the fastq files.
   | Default: ``76``

Input Directories
^^^^^^^^^^^^^^^^^

inputdir ``recommended to leave at default``
   | Base input directory.
   | Default: ``$workdir/input``

controlfolder ``recommended to leave at default``
   | Directory for all needed input files when no differential comparison. Directory for control sample input files when running differential AS event detection.
   | Default: ``$inputdir/controldir``

casefolder ``recommended to leave at default``
   | Directory only for case sample input files in case of differential AS event detection.
   | Default: ``$inputdir/casedir``

fastqdir ``recommended to leave at default``
   | Directory for fastq files.
   | Currently same as 'controlfastq'
   | Default: ``$controlfolder/fastqdir``

bamdir ``recommended to leave at default``
   | Directory for bam files.
   | Currently same as 'controlbam'
   | Default: ``$controlfolder/bamdir``

samdir ``recommended to leave at default``
   | Directory for sam files.
   | Default: ``$controlfolder/bamdir``

fastadir
   | Directory for the reference genome file
   | Default: ``$inputdir``

gtfdir
   | Directory for the annotation file file.
   | Default: ``$inputdir``

gffdir
   | Directory for gff file
   | Default: ``$inputdir``

Tool specific parameters
^^^^^^^^^^^^^^^^^^^^^^^^

bowtie_fastadir
    | Some tools require chromosome-wise fasta-inputs
    | Default: ``$inputdir/fasta_chromosomes/``



Index Parameters
^^^^^^^^^^^^^^^^

recompute_index
   | Force recompute the index even if the index with $indexname already exists.
   | Default: ``false``

indexname
   | Basename of the index (without eg. .1.bt2 for bowtie index).
   | Default: ``${fastaname}_index``

star_index
  | Folder containing a star index built with the ``$gtf`` and ``$fasta`` files (see below), used by: IRFinder, KisSplice, rMATS
  | Default: ``$workdir/index/star_index``

indexdir
   | Directory of the index.
   | Default: ``$workdir/index/${tool:-unspecific}_index``



ASimulatoR Parameters
^^^^^^^^^^^^^^^^^^^^^

asimulator_gtf
  | Name of the file in the input directory used by ASimulatoR to generate new transcripts
  | Example: ``Homo_sapiens.GRCh38.105.gtf``


Input Parameters
^^^^^^^^^^^^^^^^

fastaname
  | Name of the genome reference file (fasta format) inside ``$fastadir``.
  | Example: ``Homo_sapiens.GRCh38.dna.primary_assembly.fa``

gtfname
  | Name of annotation reference file inside ``$gffdir``.
  | Example: ``splicing_variants.gtf``

gffname
  | Name of gff reference file inside ``$gffdir``.
  | Example: ``splicing_variants.gff3``

   .. note::

    There should be no need to edit ``fasta``, ``gtf`` and ``gff`` since they just combine other parameters.

fasta
  | Full path to the reference genome file.
  | Default: ``${fastadir:-unspecific}/$fastaname``

gtf
  | Full path to the annotation file.
  | Default: ``${gtfdir:-unspecific}/$gtfname``

gff
  | Full path to the gff file.
  | Default: ``${gffdir:-unspecific}/$gffname``




Basic Mapping Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

outname ``recommended to leave at default``
    | Base name of the output files. They will usually be prefixed with the fastq file name and suffixed with ``.sam``.
    | Default: ``$tool`` (the name of the tool creating the ouput files)


.. warning::

      Something broke while changing the config file? Make sure there is no space between the variable, the equal sign and the value.

      Since these files are bash scripts, it is important to mind the syntax rules. E.g., there can't be a whitespace before and after "=".

      For example:
      | Wrong: workdir = "dockers/"
      | Right: workdir="dockers/"
