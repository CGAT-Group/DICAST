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
	
	
	echo "Running $tool unificiation..."

	tmp="${filename##*/}"
	bam_name="${tmp%%.*}"
	outdir_name="${bam_name}_output"
	echo "Looking for $tool files in $outdir/$outdir_name"

	if [[ -f "$outdir/$outdir_name/SGSeq_denovo.csv" ]]; 
	then
		unified_outdir_name="${outdir}/${outdir_name}_${tool}_unified"
		echo "Saving unified output to $unified_outdir_name"
		
		uni_tmp="/tmp/unification_tmpdir"
		mkdir $uni_tmp
		
		#Reformat SGSeq output to work with unification script
		awk -F '"' '{print $4 "\t" $6 "\t" $(NF-1)}' < ${outdir}/${outdir_name}/SGSeq_denovo.csv > $uni_tmp/SGSeq_denovo_formatted.csv

		if [ $combine_events = 0 ];
		then
			python3 /MOUNT/scripts/unified_output/output_transformer.py create --sgseq_denovo $uni_tmp/SGSeq_denovo_formatted.csv -out $unified_outdir_name -gtf $gtf
		else
			python3 /MOUNT/scripts/unified_output/output_transformer.py create --sgseq_denovo $uni_tmp/SGSeq_denovo_formatted.csv -out $unified_outdir_name -gtf $gtf -comb
		fi
	else 
		echo "Couldn't find necessary input files for unification."
	fi

	echo "Finished $tool unification for ${outdir_name}."
done
