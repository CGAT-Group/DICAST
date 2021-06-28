#!/bin/bash

tool=spladder
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

#cleaning up
trap cleaner EXIT


#make output directory
mk_outdir

#handle SAM file
handlesamfiles $differential

if [ $differential = 0 ]
then
	#list all .bam files in bamfolder with comma seperated
	while read -r bamfile; do
		bam_specific_out=$outdir/$(basename -s .bam $bamfile)_output
		mkdir -p $bam_specific_out/
		chmod -R 777 $bam_specific_out/
		echo Starting spladder in event detection mode for $(basename -s .bam $bamfile) ...
		spladder build -b $bamfile -o $bam_specific_out -a $gtf --parallel $ncores -n $read_length --output-txt-conf

		echo "Running $tool unificiation..."

		echo "Looking for $tool files in $bam_specific_out"
		unified_outdir_name="${bam_specific_out}_${tool}_unified"
		echo "Saving unified output to $unified_outdir_name"

		# Unzip spladder output files and save to new tmp directory
		uni_tmp="/tmp/unification_tmpdir"
		mkdir $uni_tmp

		for file in $bam_specific_out/*.confirmed.txt.gz ;
		do
			if [ -f "$file" ];
			then
				gunzip -c "$file" > $uni_tmp/$(basename -s .gz $file)
			fi
		done

		if [ $combine_events = 0 ];
		then
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -s $uni_tmp -out $unified_outdir_name -gtf $gtf
		else
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -s $uni_tmp -out $unified_outdir_name -gtf $gtf -comb
		fi
		echo "Finished $tool unification for ${bam_specific_out}."
	done < /tmp/controlbamlist
	
fi

if [ $differential = 1 ]
then
	caselist=$(ls -1p $casebam/*.bam | xargs echo | sed 's/ /,/g')
	controllist=$(ls -1p $controlbam/*.bam | xargs echo | sed 's/ /,/g')
	echo Starting spladder in DS mode ...
	#echo "$caselist,$controllist"
	spladder build -b "$caselist,$controllist" -o $outdir -a $gtf --parallel $ncores -n $read_length --output-txt-conf

	echo testing for differential splicing  ...
	spladder test --conditionA $caselist --conditionB $controllist -o $outdir --parallel $ncores -n $read_length --labelA "CASE" --labelB "CONTROL"
	cleaner
fi
