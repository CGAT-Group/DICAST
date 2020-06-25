#!/bin/bash

# use config and function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=sgseq


### START here ###################

#test filepath for gtf and BAM+BAM-index
test_gtf $wd/$gtf
test_bam $wd/$bamfolder/$bamfile



#make output directories
mk_outdir $tool


### Start AS event detection ###

echo compute ${tool} AS event detection...


Rscript $wd/Rscripts/SGSeq.R --gtf $wd/$gtf --path_to_bam $wd/$bamfolder/$bamfile --out $wd/$output/${tool:-unspecific}-output --cores $ncores

wait

cleaner



