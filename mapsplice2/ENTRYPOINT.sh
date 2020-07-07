#!/bin/bash

### Needed Functions

#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")
second_attempt() {
echo "paired mapping failed for ${line}. Try unpaired mapping."
python /opt/conda/bin/mapsplice.py -c /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa -x "/myvol1/index/${tool}-index/" -o /myvol1/"$tool"-output/${line##*/}${tool}/ -p 64 -1 ${line}?.fastq
}

build_index() {
#Build Genome index
# find gtf and fa files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")
mkdir -p /myvol1/index/${tool}-index
echo "Build bowtie-index..."
# Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>
#    reference_in            comma-separated list of files with ref sequences
#    ebwt_outfile_base       write Ebwt data to files with this dir/basename
# Parameters:
# --seed <int>            seed for random number generator 
bowtie-build --seed 42 /myvol1/Homo_sapiens.GRCh38.dna.primary_assembly.fa /myvol1/index/${tool}-index/
chmod -R 777 /myvol1/index/${tool}-index/
}

#cleaning up
cleaner() {
rm /myvol1/"$tool"-fastqlist
echo "script is done";exit;
}

#START here: Make a list of fastq files ############################################################################
tool=mapsplice2
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

echo "compute crac mapping..."
#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")

# python mapsplice.py [options]* -c <Reference_Sequence> -x <Bowtie_Index> -1 <Read_List1> -2 <Read_List2>
# Parameters (X means that the parameter is currently not used)
# -c <string> 	The directory containing the sequence files of reference genome.
# -x <string> 	The basename (including directory path) of Bowtie 1 index to be searched. The basename is the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / etc.
# -p / --threads <int> 	Number of threads to be used for parallel aligning. Default is 1.
# -o / --output <string> 	The directory in which MapSplice will write its output. Default is "./mapsplice_out/".
# X --gene-gtf <string> 	Gene annotation file in GTF format, used to annotate fusion junctions

python /opt/conda/bin/mapsplice.py -c /myvol1/ -x "/myvol1/index/${tool}-index/" -o /myvol1/"$tool"-output/${line##*/}${tool}/ -p 64 -1 ${line}1.fastq  -2 ${line}2.fastq

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/"$tool"-fastqlist


#make output accessible
chmod -R 777 /myvol1/"$tool"-output/

#wait for all processes to end
wait
cleaner
