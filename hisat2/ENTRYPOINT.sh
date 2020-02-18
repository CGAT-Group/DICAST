#!/bin/bash
second_attempt() {
#Attempt Unpaired alignment
hisat2 -q  -x /myvol1/hisat-output/temp/ -U "$line"?.fastq -S /myvol1/hisat-output/"${line##*/}".sam --known-splicesite-infile /myvol1/hisat-output/temp/splicesites.txt

}

#START here: Make a list of all fastq files
find /myvol1/ -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/hisat-fastqlist

#Make output directories
mkdir -p /myvol1/hisat-output/temp


#Build Genome index
cd /myvol1/hisat-output/temp/

hisat2-build /docker_main/Homo_sapiens.GRCh38.dna.primary_assembly.fa simulated_splice_data_ && python /docker_main/hisat2-2.0.0-beta/extract_splice_sites.py $(find /myvol1/ -name "*.gtf") > splicesites.txt

cd /myvol1/hisat-output/

#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
hisat2 -q  -x /myvol1/hisat-output/temp/ -1 "$line"1.fastq -2 "$line"2.fastq -S /myvol1/hisat-output/"${line##*/}".sam --known-splicesite-infile /myvol1/hisat-output/temp/splicesites.txt

#If paired end mapping fails, run unpaired mapping
trap 'second_attempt $line' ERR
done </myvol1/hisat-fastqlist

#wait for all processes to end
wait

#cleaning up
rm /myvol1/hisat-fastqlist
rm -R /myvol1/hisat-output/temp
echo "script is done"
