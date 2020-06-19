#!/bin/bash

# use config and function file
source $wd/config/asevent_config.sh
source $wd/func/asevent_func.sh
tool=sgseq


### START here ###################

#test filepath for gtf and BAM
test_gtf $wd/$gtf
test_bam $wd/$bamfolder/$sambam

#####
#readsamfiles 
#for filename in $wd/tmp/$tool-samlist do;
#	makebamfromsam $filename
#done;
#
#readbamfiles 
#####


#make output directories
mk_outdir $tool


### Start AS event detection ###

echo compute ${tool} AS event detection...

#for filename in $wd/tmp/$tool-bamlist do;
#	echo running SGSeq for $filename
#	#Rscript SGSeq.R --gtf $wd/$gtf --path_to_bam $filename --out $wd/$output/$tool-ouput --cores $ncores
#	wait
#done;

Rscript $wd/Rscripts/SGSeq.R --gtf $wd/$gtf --path_to_bam $wd/$bamfolder/$sambam --out $wd/$output/$tool-output --cores $ncores

wait

#cleaner

echo script finished.


