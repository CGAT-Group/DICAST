#!/bin/bash

tool=dSpliceType

source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

#handle SAM file
readsamfiles
for filename in $(cat $outdir/samlist)
do
	makebamfromsam $filename
done


#prepare control
echo "prepare files in $controlbam folder"
mkdir /tmp/control/
for bamfile in $(ls $controlbam/*.bam); do \
	outfile=$(basename -- ${bamfile%.*});
	echo "run bedtools for: " $bamfile;
	bedtools genomecov -ibam $bamfile > /tmp/control/${outfile}.bedgraph;
	echo "run regtools for: " $bamfile;
	regtools junctions extract $bamfile -s 0 -o /tmp/control/${outfile}.bed;
done

#prepare case
echo "prepare files in $casebam folder"
mkdir /tmp/case/
for bamfile in $(ls $casebam/*.bam); do \
	outfile=$(basename -- ${bamfile%.*});
	echo "run bedtools for: " $bamfile;
	bedtools genomecov -ibam $bamfile > /tmp/case/${outfile}.bedgraph;
	echo "run regtools for: " $bamfile;
	regtools junctions extract $bamfile -s 0 -o /tmp/case/${outfile}.bed;
done


mk_outdir

casecoverage=$(for line in $(ls -d /tmp/case/*.bedgraph); do printf "%s" ${line},; done)
controlcoverage=$(for line in $(ls -d /tmp/control/*.bedgraph); do printf "%s" ${line},; done)
casejunction=$(for line in $(ls -d /tmp/case/*.bed); do printf "%s" ${line},; done)
controljunction=$(for line in $(ls -d /tmp/control/*.bed); do printf "%s" ${line},; done)




# run analysis

#java -jar [-Xmx memory] dSpliceType.jar 
# -g <.gff> 					provide the annotation file in .gff format.
# -b1 <c1.bedgrapgh,...,cn.bedgraph>		provide read coverage (.bedgraph) files for replicates in condition 1.
# -b2 <t1.bedgraph,...,tn.bedgraph> 		provide read coverage (.bedgraph) files for replicates in condition 2.
# -j1 <c1_junc.bed,...,cn_junc.bed> 		provide junction (.bed) files for replicates in condition 1. For example, junction .bed file from Tophat or Tophat2.
# -j2 <t1_junc.bed,...,tn_junc.bed> 		provide junction (.bed) files for replicates in condition 2. For example, junction .bed file from Tophat or Tophat2.
# [options]

java -jar /docker_main/dSpliceType.jar \
	-g $gtf_dSpliceType \
	-b1 ${casecoverage%","} -b2 ${controlcoverage%","} \
	-j1 ${casejunction%","} -j2 ${controljunction%","} \
	-o $outdir

cleaner
