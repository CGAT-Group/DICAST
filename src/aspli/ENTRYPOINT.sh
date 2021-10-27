#!/bin/bash

# use config and function file
tool=aspli
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

#cleaning up
trap cleaner EXIT


#test input-files
test_gtf $gtf

#making outdir
mk_outdir
#checking for SAM files and building BAM-files with index with them
handlesamfiles $differential 2>/dev/null

echo starting ASpli...

if [ $differential = 0 ]; then
	echo Starting ASpli in AS event detection mode...
	for i in $(cat /tmp/controlbamlist)
		do j=$(basename $i)
		echo Starting ASpli run for $j -----------------------
		outdir_name=$(basename $i .bam)

		unified_outdir_name="${outdir}/${outdir_name}_${tool}_dicast_unified"
		mkdir -p /tmp/bams/$j

		ln -s $i /tmp/bams/$j/$j && ln -s $i /tmp/bams/$j/$(basename $i .bam)1.bam
		Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir/$j --bamfolder /tmp/bams/$j --differential $differential
		for k in $(find $outdir/$j/ -type f -name '*'); do mv $k $outdir/$j/$(basename $k) 2>/dev/null; rmdir $outdir/$j/*/* 2>/dev/null; rmdir $outdir/$j/* 2>/dev/null ; done
		# run unification
		echo "Running $tool unificiation..."
		anno_file="$workdir/src/ASimulatoR/out/event_annotation.tsv"
		stats_file="${unified_outdir_name}/${outdir_name}_output_${tool}_dicast_unified_comparison.txt"

		mkdir -p $unified_outdir_name
		if [ $combine_events = 0 ];
      then
        python3 /MOUNT/scripts/unified_output/output_transformer.py create --aspli_dir $outdir/$j -out $unified_outdir_name -gtf $gtf
				if [[ -f "$anno_file" ]];
				then
					echo "Running unified comparison..."
					python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0
				fi
      else
        python3 /MOUNT/scripts/unified_output/output_transformer.py create --aspli_dir $outdir/$j -out $unified_outdir_name -gtf $gtf -comb
				if [[ -f "$anno_file" ]];
				then
					echo "Running unified comparison..."
					python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0 -comb
				fi
			fi
		rm -r /tmp/bams/$j
		echo Finished ASpli run for $j -----------------------
	done
	rmdir $outdir/* 2>/dev/null
fi
if [ $differential = 1 ]; then
	echo Starting Aspli in DS analysis mode...
	Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir --casefolder $casebam --controlfolder $controlbam --differential $differential
	#Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir --bamfolder $controlbam --differential $differential
fi
