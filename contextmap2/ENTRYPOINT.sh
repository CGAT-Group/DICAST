#!/bin/bash

##VARIABLES##
fastapath=/myvol1/fasta/
indexpath=/myvol1/contextmap2-output/indices/
tool=contextmap2
export fastapath
export indexpath


#WORKING DIRECTORY
cd /myvol1

#To call contextmap as a command.
PATH="$PATH":/home/biodocker/ContextMap_Source_v2.7.9/:/home/biodocker/jre1.8.0_241/


###FUNCTIONS ###
second_attempt() {
	#Unpaired mapping command
	#...tag outputs with this flag to name it per fastqfile       "${line##*/}"
	#...address for all gtf files are                              $(find /myvol1/ -name "*.gtf")
	java -jar /home/biodocker/ContextMap_Source_v2.7.9/ContextMap_v2.7.9.jar mapper -reads ${line}?.fastq -aligner_name bowtie2 -o /myvol1/contextmap2-output -aligner_bin /home/biodocker/bin/    bowtie2 -indexer_bin /home/biodocker/bin/bowtie2-build -indices $(find ${indexpath} -name "*.bt2" | cut -f1 -d"." | tr '\n' ',' | sed s/,$//) -genome ${fastapath::-1}
}

#Build Genome index
build_index(){
	mkdir "$indexpath"
	for line in $(ls "$fastapath"); do bowtie2-build -f "$fastapath""$line" "$indexpath""$line" ; done
}


#START HERE ###########
#: Make a list of fastq files
find /myvol1/ -maxdepth 1 -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist

#make output directories
mkdir -p /myvol1/"$tool"-output/

#test filepaths for fasta and indexing
if ! test -f "${fastapath}1.fa"; then echo "check the path for the contextmap/bowtie fasta files: ${fastapath}1.fa"; exit; fi
if ! test -f "${indexpath}MT.1.bt2"; then build_index; fi

#Iterate list with paired end map command first
while read -r line; do

	echo "$line $indexpath $fastapath"
#First attempt: Paired end mapping

java -jar /home/biodocker/ContextMap_Source_v2.7.9/ContextMap_v2.7.9.jar mapper -reads ${line}1.fastq,${line}2.fastq -aligner_name bowtie2 -o /myvol1/contextmap2-output -aligner_bin /home/biodocker/bin/bowtie2 -indexer_bin /home/biodocker/bin/bowtie2-build -indices $(find ${indexpath} -name "*.bt2" | cut -f1 -d"." | tr '\n' ',' | sed s/,$//) -genome "${fastapath::-1}"

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist

#wait for all processes to end
wait

#cleaning up
rm /myvol1/"$tool"-fastqlist
rm -R /myvol1/"$tool"-output/temp 
echo "script is done"
