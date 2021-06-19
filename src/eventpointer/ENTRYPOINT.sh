#!/bin/bash

tool=eventpointer
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

##### START here ###############

#cleaning up
trap cleaner EXIT


#test filepath for gtf and BAM+BAM-index
test_gtf $gtf

#make output directory
mk_outdir
#handle SAM files
handlesamfiles $differential

#### AS event detection mode ####
if [ $differential = 0 ]; then
	echo Starting EventPointer in AS event detection mode...
	Rscript /docker_main/EventPointer.R --gtf $gtf --cores $ncores --out $outdir --bamfolder $controlbam --differential $differential
	wait
	cleaner
fi


#### DS mode ####
if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
	combined=$(combine_case_control $casebam $controlbam) #combine case and control folder into single folder
        Rscript /docker_main/EventPointer.R --gtf $gtf --cores $ncores  --out $outdir --casefolder $casefolder --controlfolder $controlfolder --differential $differential --combined $combined
	wait
	cleaner_diff
fi
