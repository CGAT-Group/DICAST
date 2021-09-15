AS detection parameters
=======================
Found in file: :guilabel:`scripts/asevent_config.sh`

This config file sets parameters that are specific to AS event detection tools only. 

Basic Parameters
^^^^^^^^^^^^^^^^

transcript
    | Fasta file for gene transcripts.
    | Default: ``$fasta`` (set in ``config.sh``, see the :ref:`documentation<input_parameters>`)

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

Differential analysis parameters (beta mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::
    For any of the parameters to be used, ``$differential`` needs to be set to ``1``.

differential
    | ``0``: Only AS event detection 
    | ``1``: All tools, which can calculate differential splicing, will do it.
    | Default: ``0``

casebam
    | Path to the case BAMs folder used for differential splicing. 
    | Default: ``$casefolder/bamdir``

casefastq
    | Path to the case fastq folder used for differential splicing.
    | Default: ``$casefolder/fastqdir``

caseprefix
    | All files in the ``$casefastq`` folder must have this prefix.
    | Example: ``sample_01``

controlbam
    | Path to the control BAMs folder used for differential splicing.
    | Default: ``$controlfolder/bamdir``

controlfastq
    | Path to the control fastq folder used for differential splicing.
    | Default: ``$controlfolder/fastqdir``

controlprefix
    | Add files in the ``controlfastq`` folder must have this prefix.
    | Example: ``sample_01``
