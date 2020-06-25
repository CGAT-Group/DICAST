#!/bin/bash

#use config & function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=eventpointer

##### START here ###############


#test filepath for gtf and BAM+BAM-index
test_gtf $wd/$gtf
test_bam $wd/$bamfolder/$bamfile


#make output directory
mk_outdir $tool


### Start AS event detection ###

echo compute $tool AS event detection

if [ $differential = 0 ]; then
        echo Starting EventPointer in AS event detection mode...
        Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $bamfolder --bamfile $bamfile --differential $differential
	wait
	cleaner
fi

if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
	combined=$(combine_case_control)
        Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores  --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $casefolder --controlfolder $controlfolder --differential $differential --combined $combined
	wait
	cleaner_diff
fi

