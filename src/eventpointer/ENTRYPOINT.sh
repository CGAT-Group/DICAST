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
	for i in $(cat /tmp/controlbamlist)
		do outdir_name=$(basename $i)
		unified_outdir_name="${outdir}/${outdir_name}_${tool}_dicast-unified"
		if [ $combine_events = 0 ];
	            then 
                python3 /MOUNT/scripts/unified_output/output_transformer.py create -e ${outdir}/${outdir_name}_output/EventsFound_RNASeq.txt -out $unified_outdir_name -gtf $gtf -comb

	            else
                python3 /MOUNT/scripts/unified_output/output_transformer.py create -e ${outdir}/${outdir_name}_output/EventsFound_RNASeq.txt -out $unified_outdir_name -gtf $gtf
	            fi
		echo Finished Eventpointer run for $outdir_name -----------------------
	done
fi


#### DS mode ####
if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
	combined=$(combine_case_control $casebam $controlbam) #combine case and control folder into single folder
        Rscript /docker_main/EventPointer.R --gtf $gtf --cores $ncores  --out $outdir --casefolder $casefolder --controlfolder $controlfolder --differential $differential --combined $combined
fi
