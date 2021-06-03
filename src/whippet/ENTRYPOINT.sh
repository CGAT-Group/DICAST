#!/bin/bash

tool=whippet
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

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

		echo "Buildiung splice graph for $filename ..."
		outdir_name=$(basename -s .bam $filename)_output
		mkdir -p $outdir/$outdir_name
		julia /docker_main/bin/whippet-index.jl --fasta $fasta --gtf $gtf --bam $filename -x $outdir/$outdir_name/graph
		echo Starting Whippet in AS event detection mode...
	        julia /docker_main/bin/whippet-quant.jl <(cat $fastqdir/*) -x $outdir/$outdir_name/graph -o $outdir/$outdir_name/whippet-out
        	cleaner

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



