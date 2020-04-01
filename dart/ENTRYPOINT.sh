#!/bin/bash
second_attempt() {
#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")
dart -i /myvol1/index/dartindex/ -f ${line}?.fastq  -o ${line}dart.sam

}

build_index() {
#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")
mkdir -p /myvol1/index/dartindex
chmod 770 /myvol1/index/dartindex
echo "compute index"
bwt_index /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa /myvol1/index/dartindex/
}

tool=dart
#START here: Make a list of fastq files
find /myvol1/ -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist
chmod 770 /myvol1/"$tool"-fastqlist

#make output directories
mkdir -p /myvol1/"$tool"-output/temp/
chmod 770 /myvol1/"$tool"-output/temp/

#test filepaths for fasta and indexing
if ! test -f "/myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa"; 
	then echo "check the path for the Homo_sapiens.GRCh* fasta files: is it under <mounted folder>/Homo_sapiens.GRCh*.fa?"; 

exit; fi
if ! test -d "myvol1/index/dartindex"; then build_index; fi




#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")
dart -i /myvol1/index/dartindex/ -f ${line}1.fastq -f2 ${line}2.fastq -o ${line}dart.sam


#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist

#wait for all processes to end
wait

#cleaning up
rm /myvol1/"$tool"-fastqlist
rm -R /myvol1/"$tool"-output/temp
echo "script is done"
