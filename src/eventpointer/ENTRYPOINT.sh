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

	# run unification
	echo "Running $tool unificiation..."
	anno_file="$workdir/src/ASimulatoR/out/event_annotation.tsv"
	for i in $(cat /tmp/controlbamlist)
		do outdir_name=$(basename ${i%%.*})
		echo "Looking for $tool files in $outdir/$outdir_name"

		if [[ -f "${outdir}/${outdir_name}_output/EventsFound_RNASeq.txt" ]];
		then
			unified_outdir_name="${outdir}/${outdir_name}_output_${tool}_dicast_unified"
			echo "Saving unified output to $unified_outdir_name"
			stats_file="${unified_outdir_name}/${outdir_name}_${tool}_dicast_unified_comparison.txt"
			mkdir -p $unified_outdir_name

			if [ $combine_events = 0 ];
      then
        python3 /MOUNT/scripts/unified_output/output_transformer.py create -e ${outdir}/${outdir_name}_output/EventsFound_RNASeq.txt -out $unified_outdir_name -gtf $gtf

				if  [[ -f "$anno_file" ]];
				then
					echo "Running unified comparison..."
					python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_output_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0
				fi
     	else
        python3 /MOUNT/scripts/unified_output/output_transformer.py create -e ${outdir}/${outdir_name}_output/EventsFound_RNASeq.txt -out $unified_outdir_name -gtf $gtf -comb
				if  [[ -f "$anno_file" ]];
				then
					echo "Running unified comparison..."
					python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_output_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0 -comb
				fi
			fi

		else
			echo "Couldn't find necessary input files for unification: ${outdir}/${outdir_name}_output/EventsFound_RNASeq.txt"
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
