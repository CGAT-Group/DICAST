#!/bin/bash

### Needed Functions

#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")
second_attempt() {
####################################################################################################################
}

#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")
build_index() {
mkdir -p /myvol1/index/dartindex
chmod 777 /myvol1/index/dartindex
echo "compute index"
# bwt_index /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa /myvol1/index/dartindex/ ###########################
}

#cleaning up
cleaner() {
rm /myvol1/"$tool"-fastqlist
echo "script is done";exit
}

#START here: Make a list of fastq files ############################################################################
tool=
find /myvol1/ -name "*fastq" -nowarn -maxdepth 1| sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist
chmod 777 /myvol1/"$tool"-fastqlist

#make output directories
mkdir -p /myvol1/index/"$tool"-index/
chmod 777 /myvol1/"$tool"-output/

#test filepaths for fasta and indexing
if ! test -e "/myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa"; 
	then echo "check the path for the Homo_sapiens.GRCh* fasta files: is it under <mounted folder>/Homo_sapiens.GRCh*.fa?";  cleaner; fi
if ! test -f "myvol1/index/${tool}-index"; then build_index; fi




#Iterate list with paired end map command first
while read -r line; do
#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile without a path     "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")#####################


#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist

#wait for all processes to end
wait
cleaner
