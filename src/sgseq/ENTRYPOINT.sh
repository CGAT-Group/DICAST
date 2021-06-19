#!/bin/bash

# use config and function file
tool=sgseq
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

### START here ###################

#cleaning up
trap cleaner EXIT

#tests
test_gtf $gtf


#handle sam files
handlesamfiles 0


#build bamlist
readbamfiles $controlbam controlbamlist
#Make output directories
mk_outdir
### Start AS event detection ###

echo compute ${tool} AS event detection...

for filename in $(cat /tmp/controlbamlist)
do
	echo Starting SGSeq for $filename ...
	sample_out=$(mk_sample_out $filename)
	Rscript /docker_main/SGSeq.R --gtf $gtf --path_to_bam $filename --out $sample_out --cores $ncores
	wait
done
