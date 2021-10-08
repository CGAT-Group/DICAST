Mapping parameters
==================
Found in file: :guilabel:`scripts/mapping_config.sh`

This config file sets parameters that are specific to mapping tools only. 

Basic Parameters
^^^^^^^^^^^^^^^^

outname
    | Base name of the output files. They will usually be prefixed with the fastq file name and suffixed with ``.sam``.
    | Default: ``$tool`` (the name of the tool creating the ouput files)

Index Parameters
^^^^^^^^^^^^^^^^

recompute_index
    | Force recompute the index even if the index with $indexname already exists.
    | Default: ``false``

indexname
    | Basename of the index (without eg. .1.bt2 for bowtie index).
    | Default: ``${fastaname}_index``

indexdir
    | Directory of the index.
    | Default: ``$workdir/index/${tool:-unspecific}_index``

Tool specific parameters
^^^^^^^^^^^^^^^^^^^^^^^^

bowtie_fastadir
    | Some tools require chromosome-wise fasta-inputs
    | Default: ``$inputdir/fasta_chromosomes/``
