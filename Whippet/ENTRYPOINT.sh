#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=whippet


mk_outdir $tool
test_gtf $wd/$gtf
test_fasta $wd/$fasta
out=$wd/$output/${tool:-unspecific}-output #local variable for output folder


readsamfiles
for filename in $(cat $wd/$output/${tool:-unspecific}-output/samlist)
do
        makebamfromsam $filename
done


echo "Merging and indexing BAM files of $wd/$bamfolder to use for splice graph building..."
readbamfiles
bamlist=$(cat $out/bamlist | xargs echo)
mkdir -p $out/tmp
samtools merge $out/tmp/merged_bam.bam $bamlist
samtools sort $out/tmp/merged_bam.bam $out/tmp/merged_bam_sorted
samtools rmdup -S $out/tmp/merged_bam_sorted.bam $out/tmp/merged_bam_sorted_rmdup.bam
samtools index $out/tmp/merged_bam_sorted_rmdup.bam

echo ----------Starting Whippet---------

echo Building splice graph...
julia /docker_main/bin/whippet-index.jl --fasta $wd/$fasta --gtf $wd/$gtf --bam $out/tmp/merged_bam_sorted_rmdup.bam -x $out/graph 

if [ $differential = 0  ]; then
	echo Starting Whippet in AS event detection mode...
	julia /docker_main/bin/whippet-quant.jl <(cat $wd/$fastqfolder/*) -x $out/graph -o $out/whippet-out
	cleaner
fi

if [ $differential = 1  ]; then
	echo Starting Whippet in differential analysis mode...
	julia /docker_main/bin/whippet-quant.jl <(cat $wd/$casefastq/*) -x $out/graph -o $out/whippet-case-out
	julia /docker_main/bin/whippet-quant.jl <(cat $wd/$controlfastq/*) -x $out/graph -o $out/whippet-control-out

	echo "Comparing case & control output..."
	julia /docker_main/bin/whippet-delta.jl -a $out/whippet-case-out.psi.gz, -b $out/whippet-control-out.psi.gz, -o $out/differential-analysis-out
	cleaner
fi
