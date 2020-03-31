#!/bin/bash
second_attempt() {
#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")



}
tool=
#START here: Make a list of fastq files
find /myvol1/ -name "*fastq" -nowarn -maxdepth 1| sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist

#make output directories
mkdir -p /myvol1/"$tool"-output/temp/

#test filepaths for fasta and indexing
#if ! test -f "/myvol1/Homo_sapiens.GRCh*.fa"; then echo "check the path for the Homo_sapiens.GRCh* fasta files: is it under <mounted folder>/Homo_sapiens.GRCh*.fa?"; exit; fi
#if ! test -f "myvol1/star-output/temp/genomedir/SAindex"; then build_index; fi


#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")


#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")


#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist

#wait for all processes to end
wait

#cleaning up
rm /myvol1/"$tool"-fastqlist
rm -R /myvol1/"$tool"-output/temp
echo "script is done"
