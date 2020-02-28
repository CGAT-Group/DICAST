#!/bin/bash

cd /myvol1

PATH="$PATH":/home/biodocker/ContextMap_Source_v2.7.9/:/home/biodocker/jre1.8.0_241/

#/** second_attempt() {
#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")


#}


tool=contextmap2
#START here: Make a list of fastq files
find /myvol1/ -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist

#make output directories
mkdir -p /myvol1/"$tool"-output/temp/

#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*Homo*.fa")  $(find /myvol1/ -name "*.gtf")
java -jar /home/biodocker/ContextMap_Source_v2.7.9/ContextMap_v2.7.9.jar indexer -fasta /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa -o /myvol1/contextmap2-output/temp/

#Iterate list with paired end map command first
while read -r line; do

	echo "$line"
#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")

#If paired end mapping fails, run unpaired mapping.
#trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist

#wait for all processes to end
wait

#cleaning up
rm /myvol1/"$tool"-fastqlist
rm -R /myvol1/"$tool"-output/temp
echo "script is done"

