#!/bin/bash

### Needed Functions

#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")
second_attempt() {
echo "paired mapping failed for ${line}. Try unpaired mapping."
gsnap --db "Homo_sapiens.GRCh38" --dir "/myvol1/index/${tool}-index" --output-file /myvol1/"$tool"-output/${line##*/}${tool}.sam --format sam --nthreads 64 ${line}?.fastq

}

build_index() {
#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")
mkdir -p /myvol1/index/${tool}-index
echo "compute gmap database"

# Usage: gmap_build [options...] -d <genomename> <fasta_files>
# -D, --dir=STRING          Destination directory for installation (defaults to gmapdb directory specified at configure time)
# -d, --db=STRING           Genome name

gmap_build -d "Homo_sapiens.GRCh38" -dir "/myvol1/index/${tool}-index" $(ls /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa)

chmod -R 777 /myvol1/index/${tool}-index/
}

#cleaning up
cleaner() {
rm /myvol1/"$tool"-fastqlist
echo "script is done";exit;
}

#START here: Make a list of fastq files ############################################################################
tool=gsnap
find /myvol1/ -name "*fastq" -nowarn -maxdepth 1| sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/"$tool"-fastqlist
chmod 777 /myvol1/"$tool"-fastqlist

#test filepaths for fasta and indexing
if ! test -f "/myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa"; 
	then echo "check the path for the Homo_sapiens.GRCh* fasta files: is it under <mounted folder>/Homo_sapiens.GRCh*.fa?";
	cleaner;
	exit; fi
#test if index folder exists
if ! test -d "/myvol1/index/${tool}-index"; then mkdir -p /myvol1/index/${tool}-index; fi
# test if index is already computed
if ! [ "$(ls -a /myvol1/index/${tool}-index)" ]; then build_index; fi

#make output directories
mkdir -p /myvol1/"$tool"-output/
chmod -R 777 /myvol1/"$tool"-output/

echo "compute $tool mapping..."
#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")

# Usage: gsnap [OPTIONS...] <FASTA file>, or
# Parameters (X means that the parameter is currently not used)
# --db	Genome database
# --dir Genome directory
# --nthreads	Number of worker threads
# --output-file File name for a single stream of output results.
# --format=STRING Another format type, other than default. Currently implemented: sam, m8 (BLAST tabular format)
# X --ordered      Print output in same order as input (relevant only if there is more than one worker thread)
# ...

gsnap --db "Homo_sapiens.GRCh38" --dir "/myvol1/index/${tool}-index" --output-file /myvol1/"$tool"-output/${line##*/}${tool}.sam --format sam --nthreads 64 ${line}1.fastq ${line}2.fastq

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist

#make output accessible
chmod -R 777 /myvol1/"$tool"-output/

#wait for all processes to end
wait
cleaner
