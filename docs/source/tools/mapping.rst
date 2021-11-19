Mapping tools
=============

.. note::

  Most mapping tools need an index file of the reference genome for mapping. The computation of these index files can take a long time.
  Our mapping tool scripts check if there already is an index for the respective tool and only build it, if it is not found.
  If you face any index related errors, either set the parameter ``$recompute_index=True`` or delete the old index to recalculate it.

Mapping Input Files
^^^^^^^^^^^^^^^^^^^^^^

.. tip::
  The paths assume you are using our suggested :doc:`input structure <../get-started/setup/input>`.
  Example input files you can find in our :doc:`examples section<../get-started/run/example-files>`.

You can find the required input files in the tool-specific documentation.

.. _fastqMapping:

fastq
  Fastq files for paired end mapping. The directories are separated in ``controldir`` and ``casedir``. The controldir is the default folder for all analyses. The casedir is only used for differential splicing analysis.

  .. code-block:: bash

    input/fastq/controldir/*yourFastqFile1*_1.fastq
    input/fastq/controldir/*yourFastqFile1*_2.fastq
    input/fastq/controldir/*yourFastqFile2*_1.fastq
    input/fastq/controldir/*yourFastqFile2*_2.fastq
    . . .

.. _fastaMapping:

fasta
  The fasta reference for your organism. Mapping tools usually only need it for indexing (see tool specific documentation).

  .. code-block:: bash

    input/*yourFastaFile*.fa

.. _gtfMapping:

gtf
  annotation reference file.

  .. code-block:: bash

    input/*yourGTFfile*.gtf

.. _bowtie_fastadirMapping:

bowtie_fastadir
  Only needed by some tools. Chromosome-wise fasta files for your organism to build an index with bowtie.

  .. code-block:: bash

    input/bowtie_fastadir/1.fa
    input/bowtie_fastadir/2.fa
    input/bowtie_fastadir/3.fa
    input/bowtie_fastadir/4.fa
    . . .
    input/bowtie_fastadir/X.fa
    input/bowtie_fastadir/Y.fa

Optional: Index
  Tool specific index file(s). If no index file is found in the index folder it will be built the first time you run the tool.
  **This might take some time.** If you want to provide your own index please make sure it is in the correct format and file names. Since the index is usually built based on the fasta reference we recommend to name the index based on the fasta reference (default). You can change the ``indexname`` variable in the config script.

  .. code-block:: bash

    index/*toolname*-index/*yourIndexBaseName*


Parameters
^^^^^^^^^^^^^

To provide a fair baseline while maintaining easy usability, per default we run the tools with their default variables. The default parameters can be changed by editing the ENTRYPOINT.sh scripts of each tool. The variables used by mapping ENTRYPOINT.sh scripts can be set in the ``config.sh`` and ``mapping_config.sh`` files in the ``scripts`` folder. For a usual analysis you should not need to change these parameters.


.. toctree::
   :maxdepth: 2
   :hidden:

   mapping/bbmap
   mapping/contextmap
   mapping/crac
   mapping/dart
   mapping/gsnap
   mapping/hisat
   mapping/mapsplice
   mapping/minimap
   mapping/segemehl
   mapping/star
   mapping/subjunc
