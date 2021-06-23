### Set parameters and add them to the params list ----
### For more information, see https://github.com/biomedbigdata/ASimulatoR

### Don't change input and output unless you know what you are doing!
args = commandArgs(trailingOnly = TRUE)
input = args[1]
output = args[2]
###

ncores = 20
multi_events_per_exon = T
probs_as_freq = F
error_rate = 0.001
readlen = 76
max_genes = 100 # NULL 
seq_depth = 10000 # 2e08
num_reps = c(1,1)
as_events = c('es', 'mes', 'ir', 'a3', 'a5', 'afe', 'ale', 'mee')

as_combs = combn(as_events, 2, FUN = function(...) paste(..., collapse = ','))
event_probs = rep(1/(length(as_combs) + 1), length(as_combs))
names(event_probs) = as_combs

outdir = sprintf(
    'run_result',
    output,
	gsub('[ ]+?', '-', Sys.time()),
    ifelse(is.null(max_genes), 0, max_genes),
    seq_depth,
    error_rate,
    readlen,
    multi_events_per_exon,
    probs_as_freq
  )

params = list(
  ncores = ncores,
  input_dir = input,
  event_probs = event_probs,
  outdir = outdir,
  seq_depth = seq_depth,
  max_genes = max_genes,
  error_rate = error_rate,
  readlen = readlen,
  multi_events_per_exon = multi_events_per_exon,
  probs_as_freq = probs_as_freq,
  num_reps = num_reps
)


### run simulator ----
library(ASimulatoR)
do.call(simulate_alternative_splicing, params)


### copy this script to output folder for reproducibility ----
rscript = sub('--file=', '', commandArgs()[4], fixed = T)
bn_rscript = basename(rscript)
file.copy(rscript, file.path(outdir, bn_rscript))
