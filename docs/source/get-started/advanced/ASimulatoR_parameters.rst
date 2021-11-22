ASimulatoR config
===================

Found in file: :guilabel:`scripts/ASimulatoR_config.R`

Parameters are also explained in the following git `github/biomedbigdata/ASimulatoR <https://github.com/biomedbigdata/ASimulatoR>`_

ncores
  | Number of cores used by ASimulatoR
  | Within dicast, the max number of cores supplied to Snakemake, is as much is used within DICAST's pipeline. ( See :doc:`Snakemake parameters <snakemake_parameters>`)

multi_events_per_exon
  | ``T`` or ``F``
  | Should each exon be treated as a target for only one Alternative Splicing event or would you like to see events like Multiple Exon Skipping events, Alternative Last/First Exon + Exon Skipping events?

probs_as_freq
  | ``T`` or ``F``
  | Default:
  | ``F``: if probs_as_freq was FALSE, a random number would be drawn for each event-superset combination and only if it was smaller than 1/9 the AS event would be created
  | ``T``: The exon supersets are partitioned corresponding to the event_prob parameter.

error_rate
  | Default: 0.001
  | In the uniform error model, probability that the sequencer records the wrong nucleotide at any given base.

readlen
  | Read Length
  | Default is 76

max_genes
  | define the number of genes you want to work with. If you want all exons, do not specify this parameter or set it to NULL

seq_depth
  | Sequencing depth of simulated experiment

num_reps
  | define, how many groups and samples per group you analyze. Here we create a small experiment with two groups with one sample per group:

as_events
  | make a list in R with the following set or a subset of the following:
  | c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee')

as_combs
  | Combinations of AS events desired in the simulated dataset.

event_probs
  | Event probabilities of AS events within the simulated dataset
