#!/bin/bash
second_attempt() {
#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")
echo "paired mapping failed. Try unpaired mapping."
dart -i /myvol1/index/dartindex/ -f ${line}?.fastq  -o /myvol1/"$tool"-output/${line##*/}dart.sam
}

build_index() {
#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")
mkdir -p /myvol1/index/${tool}-index
echo "compute index"
bwt_index /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa /myvol1/index/${tool}-index/
chmod -R 777 /myvol1/index/${tool}-index
}

cleaner() {
rm /myvol1/"$tool"-fastqlist
echo "script is done";exit;
}

tool=dart
#START here: Make a list of fastq files
find /myvol1/ -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist
chmod -R 777 /myvol1/"$tool"-fastqlist



#test filepaths for fasta and indexing
if ! test -e "/myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa"; 
	then echo "check the path for the Homo_sapiens.GRCh* fasta files: is it under <mounted folder>/Homo_sapiens.GRCh*.fa?"; 
	cleaner;
	exit; fi
if ! test -d "/myvol1/index/dartindex"; then build_index; fi

#make output directories
mkdir -p /myvol1/"$tool"-output/
chmod -R 777 /myvol1/"$tool"-output/


#Iterate list with paired end map command first
while read -r line; do
#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")
dart -i /myvol1/index/dartindex/ -f ${line}1.fastq -f2 ${line}2.fastq -o /myvol1/"$tool"-output/${line##*/}dart.sam
#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist


#make output accessible
chmod -R 777 /myvol1/"$tool"-output/

#wait for all processes to end
wait
cleaner

