Configuring DICAST
====================

Before running DICAST, please take some time to configure it for your first run.  

DICAST is best :doc:`run with the GUI <run/gui>`, which automates the configuration of :guilabel:`scripts/snakemake/snakemake_config.yaml`, :guilabel:`scripts/config.sh` & :guilabel:`scripts/asevent_config.sh`,
for a quick run of DICAST without simulated data, on your experiments. 

DICAST can be run :doc:`via CLI <run/dicast-cli>`, however, this feature is currently in development.

If you'd like to modify the **Simulated dataset**, please modify :guilabel:`scripts/ASimulatoR_config.R` (See :doc:`ASimulatoR Parameters <config/ASimulatoR_parameters>`)

If you'd like to run DICAST with just one tool :doc:`via docker <run/one_tool>`, then the files you need to modify are: :guilabel:`scripts/config.sh`, :guilabel:`scripts/asevent_config.sh`

It's recommended to take a closer look at the config files on disk before your first run.

The following files are all the configuration files found in DICAST:
  | :guilabel:`scripts/snakemake/snakemake_config.yaml`,
  | :guilabel:`scripts/ASimulatoR_config.R`,
  | :guilabel:`scripts/config.sh` &
  | :guilabel:`scripts/asevent_config.sh`.


Read the full reference here
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  .. toctree::
    :maxdepth: 2

    config/snakemake_parameters
    config/ASimulatoR_parameters
    config/general_parameters
    config/as_parameters


Most frequent changes to configurations:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

	The following explanations assume that you use our directory structure as described in :doc:`Directory Structure <setup/input>`.

1. :guilabel:`scripts/snakemake/snakemake_config.yaml`
*********************************************************

The following are the snakemake parameters that you're most likely to change for a :doc:`CLI run: <run/dicast-cli>`.

Possible_overwrite_acknowledge:
  | do: false
  | change to true. This is set to false after every run to prevent overwriting of output files

Mapping_tools:
  |    What_tools_to_run: '<insert name of mapping tools to run, separated by spaces>'
  | pick one of the following :guilabel:`bbmap` :guilabel:`contextmap` :guilabel:`crac` :guilabel:`dart` :guilabel:`gsnap` :guilabel:`hisat` :guilabel:`mapsplice` :guilabel:`minimap` :guilabel:`segemehl` :guilabel:`star` :guilabel:`subjunc`
  | Example: to some two tools: 'minimap star'
  | Example: to run all tools: 'bbmap contextmap crac dart gsnap hisat mapsplice minimap segemehl star subjunc'
  | Example: to run one tool: 'star'

Alternative_splicing_detection_tools:
  |    What_tools_to_run: '<insert name of Alternative Splicing tools to run, separated by spaces>'
  |  pick one of the following :guilabel:`asgal` :guilabel:`aspli` :guilabel:`eventpointer` :guilabel:`irfinder` :guilabel:`majiq` :guilabel:`sgseq` :guilabel:`spladder` :guilabel:`whippet`
  | Example: to run all tools: 'asgal aspli eventpointer irfinder majiq sgseq spladder whippet'
  | Example: to some two tools: 'eventpointer whippet'
  | Example: to run one tool: 'whippet'


2. :guilabel:`scripts/config.sh`
*********************************

The following are basic parameters that you are most likely to change on the GUI and in the file :doc:`scipts/config.sh <config/general_parameters>`.


.. warning::

  Since these files are bash scripts, it is important to mind the syntax rules. E.g., there can't be a whitespace before and after "=".


Basic Parameters
^^^^^^^^^^^^^^^^

ncores
	| if you want to use more cores for each tool within snakemake. (not the same as total cores available for snakemake ``-j 2``)


Input Parameters
^^^^^^^^^^^^^^^^
asimulator_gtf
	| the genome gtf annotation that you use to simulate the data. Default: 'Homo_sapiens.GRCh38.104.gtf'.
fastaname
	| the genome reference file. Default: 'Homo_sapiens.GRCh38.dna.primary_assembly.fa'.
gtfname
	| the genome gtf annotation that you use for mapping and alternative splicing analysis. If you’re using ASimulatoR, leave this as ASimulatoR.gtf.
gffname
	| the genome gff3 annotation that you use for mapping and alternative splicing analysis. If you’re using ASimulatoR, leave this as ASimulatoR.gff3.

The reference genome, annotation file and gff3 files could be downloaded from `Ensembl <http://ftp.ensembl.org/pub/release-104/>`_.
