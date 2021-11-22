AS detection parameters
=======================
Found in file: :guilabel:`scripts/asevent_config.sh`

This config file sets parameters that are specific to AS event detection tools only.

Basic Parameters
^^^^^^^^^^^^^^^^

transcript
    | Fasta file for gene transcripts.
    | Default: ``$fasta`` (set in ``config.sh``, see the :doc:`General Parameters <general_parameters>`)

star_alignment_files
    | Path to the folder containing star alignment files (\*.SJ.)
    | Default: ``$workdir/output/star-output``

.. note::
    | We support only paired RNA-Seq - fastq files have to be in pairs.
    | Set the suffixes parameters (including the file extension) for all fastq pairs (e.g. ``_1.fastq`` and ``_2.fastq``).

fastqpair1suffix
    | Suffix for the first file of the fastq pair.
    | Example: ``_1.fastq``

fastqpair2suffix
    | Suffix for the second file of the fastq pair.
    | Example: ``_2.fastq``

use_bam_input_files
    | Determines what kind of input to use: ``1`` for bam files, ``0`` for fastq files.
    | Default: ``0``

combine_events
    | Events such as Multiple Exon Skipping should be represented as such, instead of individual exon skipping events.
    | Default: ``1``
