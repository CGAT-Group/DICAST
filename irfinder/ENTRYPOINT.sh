#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=irfinder


mk_outdir $tool
out=${wd}/${output}/${tool:-unspecific}-output  #local variable for readability

#######################################
### build IRFinder custom reference ###
#######################################

#test for transcript.gtf and genome.fa
gtf=$wd/$gtf
fasta=$wd/$fasta
test_gtf $gtf
test_fasta $fasta


#move both files into tmp folder
echo linking annotation files into reference folder...
mkdir -p $out/index
ln -sf $gtf $out/index/transcripts.gtf
ln -sf $fasta $out/index/genome.fa


#build reference
echo building reference...
IRFinder -m BuildRefProcess -r $out/index 
wait
echo reference built, moving on...


#####################
#### quantify IR ####
#####################

#readfastqs $fastqfolder
#fastqs=$(cat $out/fastqlist)

echo looking for intron retention events...
readbamfiles
bams=$(cat $out/bamlist)
for bam in $bams
do
       	IRFinder -m BAM -r $out/index -d $out $bam
       	wait
done

 
cleaner
