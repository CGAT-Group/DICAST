Snakemake parameters
=====================

NEEDS edit
Found in file: :guilabel:`scripts/snakemake/snakemake_config.yaml`

These parameters are either set in the GUI, or if you're running DICAST via cli, these parameters determine properties your DICAST run.

Possible_overwrite_acknowledge:
  | When running DICAST first, an :guilabel:`output` directory is created. When your run is finished, please rename the output directory, if you want to save this output; so that you don't overwrite outputs with the second run of DICAST.
  | ASimulatoR files such as src/ASimulatoR/out/event_annotation.tsv are also overwritten between runs, if ASimulatoR is run again.
  | :guilabel:`true` or :guilabel:`false`
  | :guilabel:`true`: DICAST runs uninterrupted.
  | :guilabel:`false`: DICAST run is interrupted until `true`

ASimulatoR:
  | do: :guilabel:`True` / :guilabel:`False`
  | Run ASimulatoR with the configs as stored in file: ::guilabel:`scripts/snakemake/snakemake_config`
  | ( See :doc:`ASimulatoR config <ASimulatoR_parameters>`)

Mapping_tools:
  |    What_tools_to_run: '<insert name of mapping tools to run, separated by spaces>'
  | pick one of the following :guilabel:`bbmap` :guilabel:`contextmap` :guilabel:`crac` :guilabel:`dart` :guilabel:`gsnap` :guilabel:`hisat` :guilabel:`mapsplice` :guilabel:`minimap` :guilabel:`segemehl` :guilabel:`star` :guilabel:`subjunc`
  | Example: to run all tools: 'bbmap contextmap crac dart gsnap hisat mapsplice minimap segemehl star subjunc'
  | Example: to some two tools: 'minimap star'
  | Example: to run one tool: 'star'

Alternative_splicing_detection_tools:
  |    What_tools_to_run: '<insert name of Alternative Splicing tools to run, separated by spaces>'
  |  pick one of the following :guilabel:`asgal` :guilabel:`aspli` :guilabel:`eventpointer` :guilabel:`irfinder` :guilabel:`majiq` :guilabel:`sgseq` :guilabel:`spladder` :guilabel:`whippet`
  | Example: to run all tools: 'asgal aspli eventpointer irfinder majiq sgseq spladder whippet'
  | Example: to some two tools: 'eventpointer whippet'
  | Example: to run one tool: 'whippet'
