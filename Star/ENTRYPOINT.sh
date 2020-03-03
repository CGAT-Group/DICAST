#!/bin/bash
second_attempt() {
#Unpaired mapping command
/docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir /myvol1/star-output/temp/genomedir --outFileNamePrefix /myvol1/star-output/Star_mapped_"${line##*/}" --sjdbGTFfile $(find /myvol1/ -maxdepth 1 -name "*.gtf")  --twopassMode Basic --runThreadN 60 --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn "$line"?.fastq

}
#START here: Make a list of fastq files

find /myvol1/ -maxdepth 1 -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/star-fastqlist

#make output directories
mkdir -p /myvol1/star-output/temp/genomedir
cd /myvol1/star-output

#Build Genome index
/docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /myvol1/star-output/temp/genomedir --genomeFastaFiles $(find /myvol1/  -maxdepth 1 -name "*.fa") --runThreadN 60 --sjdbGTFfile $(find /myvol1/  -maxdepth 1 -name "*.gtf") --sjdbOverhang 100


#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
/docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir /myvol1/star-output/temp/genomedir --outFileNamePrefix /myvol1/star-output/Star_mapped_"${line##*/}" --sjdbGTFfile $(find /myvol1/ -maxdepth 1  -name "*.gtf")  --twopassMode Basic --runThreadN 60 --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn "$line"1.fastq "$line"2.fastq

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/star-fastqlist

#wait for all processes to end
wait

#cleaning up
#rm /myvol1/star-fastqlist
#rm -R /myvol1/star-output/temp
echo "script is done"
