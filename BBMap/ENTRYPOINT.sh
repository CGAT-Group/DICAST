#!/bin/bash

### Tool-specific variables ###
tool=BBMap

# use confic and function file
source /myvol1/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	/bbmap/bbmap.sh \
		in=${line}?.fastq \
		outm=$out/${line##*/}${tool}_mapped.sam \
		outu=$out/${line##*/}${tool}_unmapped.sam
}

build_index() {
	mkdir -p $indexdir/$index
	echo "compute index ..."
	/bbmap/bbmap.sh ref=$(ls $inputdir/$fasta) path=$indexdir/$index overwrite=false
	chmod -R 777 $indexdir
	echo "Index is now saved at $indexdir/$index"
}

### START here ############################################################################

# test filepaths
test_fasta

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -d $indexdir/$index; then build_index; fi fi

#make list of fastq files
mk_fastqlist

#make output directories
mk_outdir

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
	#First attempt: Paired end mapping
	#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
	#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")

	/bbmap/bbmap.sh \
		in=${line}1.fastq \
		in2=${line}2.fastq \
		ref=$(ls $inputdir/$fasta) \
		path=$indexdir/$index \
		outm=$out/${line##*/}${tool}_mapped.sam \
		outu=$out/${line##*/}${tool}_unmapped.sam

	#If paired end mapping fails, run unpaired mapping.
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner

