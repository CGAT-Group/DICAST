#!/bin/bash

tool=whippet
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

#cleaning up
trap cleaner EXIT

#make output directory
mk_outdir
#test reference files
test_gtf $gtf
test_fasta $fasta

#handle sam files
handlesamfiles $differential

#echo "Merging and indexing BAM files of $controldir to use for splice graph building..."
bamlist=$(cat /tmp/casebamlist /tmp/controlbamlist | xargs echo)
mkdir -p $outdir/tmp

if [ $differential = 0  ]; then
	for filename in $(cat /tmp/controlbamlist)
	do
		# no need to index each BAM file --> index should already exist!
		#samtools sort $filename -o $outdir/tmp/$filename.sorted.bam
		#samtools rmdup -S $outdir/tmp/$filename.sorted.bam $outdir/tmp/$filename.sorted.rmdup.bam
		#samtools index $outdir/tmp/$filename.sorted.rmdup.bam

		echo "Building splice graph for $filename ..."
		outdir_name=$(basename -s .bam $filename)_output
		mkdir -p $outdir/$outdir_name
		julia /docker_main/bin/whippet-index.jl --fasta $fasta --gtf $gtf --bam $filename -x $outdir/$outdir_name/graph
		echo Starting Whippet in AS event detection mode...
		julia /docker_main/bin/whippet-quant.jl <(cat $fastqdir/*) -x $outdir/$outdir_name/graph -o $outdir/$outdir_name/whippet-out
		cleaner

		echo "Running $tool unificiation..."

		echo "Looking for $tool files in $outdir/$outdir_name"
		unified_outdir_name="${outdir}/${outdir_name}_${tool}_dicast_unify"
		mkdir -p $unified_outdir_name
		echo "Saving unified output to $unified_outdir_name"

		# Unzip whippet output files and save to new tmp directory
		uni_tmp="/tmp/unification_tmpdir"
		mkdir -p $uni_tmp
		gunzip -c $outdir/$outdir_name/whippet-out.psi.gz > $uni_tmp/whippet-out.psi

		anno_file="$workdir/src/ASimulatoR/out/event_annotation.tsv"
		stats_file="${unified_outdir_name}/${outdir_name}_${tool}_dicast_unify_comparison.txt"

		if [ $combine_events = 0 ];
		then
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -w $uni_tmp/whippet-out.psi -out $unified_outdir_name -gtf $gtf
			if [[ -f "$anno_file" ]];
			then
				echo "Running unified comparison..."
				python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_${tool}_dicast_unify.out -gtf $gtf -stats $stats_file -s -t 0
			fi
		else
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -w $uni_tmp/whippet-out.psi -out $unified_outdir_name -gtf $gtf -comb

			if [[ -f "$anno_file" ]];
			then
				echo "Running unified comparison..."
				python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_${tool}_dicast_unify.out -gtf $gtf -stats $stats_file -s -t 0 -comb
			fi
		fi
		echo "Finished $tool unification for ${outdir_name}."
	done
fi

if [ $differential = 1  ]; then
	echo "Merging and indexing BAM files of $controldir to use for splice graph building..."
	samtools merge $outdir/tmp/merged_bam.bam $bamlist
	samtools sort $outdir/tmp/merged_bam.bam -o $outdir/tmp/merged_bam_sorted.bam
	samtools rmdup -S $outdir/tmp/merged_bam_sorted.bam $outdir/tmp/merged_bam_sorted_rmdup.bam
	samtools index $outdir/tmp/merged_bam_sorted_rmdup.bam

	echo Building splice graph for merged bam file...
	julia /docker_main/bin/whippet-index.jl --fasta $fasta --gtf $gtf --bam $outdir/tmp/merged_bam_sorted_rmdup.bam -x $outdir/graph

        echo Starting Whippet in differential analysis mode...
        julia /docker_main/bin/whippet-quant.jl <(cat $casefastq/*) -x $outdir/graph -o $outdir/whippet-case-out
        julia /docker_main/bin/whippet-quant.jl <(cat $controlfastq/*) -x $outdir/graph -o $outdir/whippet-control-out

        echo "Comparing case & control output..."
        julia /docker_main/bin/whippet-delta.jl -a $outdir/whippet-case-out.psi.gz, -b $outdir/whippet-control-out.psi.gz, -o $outdir/differential-analysis-out
        cleaner
fi
