#!/bin/bash

tool=dSpliceType

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh

# create bedgraph files (read coverage)
bedtools genomecov -ibam [condition_1.bam] > [condition_1.bedgraph] #for each sample

# create bed files (junctions)
regtools junctions extract [condition_1.bam] -s 0 > [condition_1.bed] #for each sample for each condition separately

# run analysis

#java -jar [-Xmx memory] dSpliceType.jar 
# -g <.gff> 					provide the annotation file in .gff format.
# -b1 <c1.bedgrapgh,...,cn.bedgraph>		provide read coverage (.bedgraph) files for replicates in condition 1.
# -b2 <t1.bedgraph,...,tn.bedgraph> 		provide read coverage (.bedgraph) files for replicates in condition 2.
# -j1 <c1_junc.bed,...,cn_junc.bed> 		provide junction (.bed) files for replicates in condition 1. For example, junction .bed file from Tophat or Tophat2.
# -j2 <t1_junc.bed,...,tn_junc.bed> [options]	provide junction (.bed) files for replicates in condition 2. For example, junction .bed file from Tophat or Tophat2.

java -jar /opt/dSpliceType.jar \
	-g /myvol1/$gff \
	-b1 \
	-b2 \
	-j1 \
	-j2 \
