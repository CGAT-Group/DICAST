#!/bin/bash

tool=whippet
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh


#mk_outdir $tool
test_gtf $gtf
test_fasta $fasta

#handle sam files
handlesamfiles $differential

echo "Merging and indexing BAM files of $bamdir to use for splice graph building..."
readbamfiles
bamlist=$(cat $outdir/bamlist | xargs echo)
mkdir -p $outdir/tmp

samtools merge $outdir/tmp/merged_bam.bam $bamlist
samtools sort $outdir/tmp/merged_bam.bam -o $outdir/tmp/merged_bam_sorted.bam
samtools rmdup -S $outdir/tmp/merged_bam_sorted.bam $outdir/tmp/merged_bam_sorted_rmdup.bam
samtools index $outdir/tmp/merged_bam_sorted_rmdup.bam

echo ----------Starting Whippet---------

echo Building splice graph...
julia /docker_main/bin/whippet-index.jl --fasta $fasta --gtf $gtf --bam $outdir/tmp/merged_bam_sorted_rmdup.bam -x $outdir/graph 

if [ $differential = 0  ]; then
	echo Starting Whippet in AS event detection mode...
	julia /docker_main/bin/whippet-quant.jl <(cat $fastqdir/*) -x $outdir/graph -o $outdir/whippet-out
	cleaner
fi

if [ $differential = 1  ]; then
	echo Starting Whippet in differential analysis mode...
	julia /docker_main/bin/whippet-quant.jl <(cat $casefastq/*) -x $outdir/graph -o $outdir/whippet-case-out
	julia /docker_main/bin/whippet-quant.jl <(cat $controlfastq/*) -x $outdir/graph -o $outdir/whippet-control-out

	echo "Comparing case & control output..."
	julia /docker_main/bin/whippet-delta.jl -a $outdir/whippet-case-out.psi.gz, -b $outdir/whippet-control-out.psi.gz, -o $outdir/differential-analysis-out
	cleaner
fi
