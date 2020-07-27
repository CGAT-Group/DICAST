#!/bin/bash

### Tool-specific variables ###
tool=hisat

# use confic and function file
source /myvol1/config/mapping_config.sh
source /myvol1/func/mapping_func.sh

#Update PATH
PATH=$PATH:/docker_main/hisat2-2.0.0-beta

### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	hisat2 -q \
	-x /$wd/index/$tool-index/$index \
	-U "$line"?.fastq \
	-S $out/${line##*/}${tool}.sam \
	--known-splicesite-infile /$wd/index/$tool-index/$index.splicesites.txt
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p /$wd/index/$tool-index
	echo "compute index ..."
	hisat2-build  $(ls $inputdir/$fasta) /$wd/index/$tool-index/$index && python /docker_main/hisat2-2.0.0-beta/extract_splice_sites.py $(ls /$wd/$inputdir/$gtf) > /$wd/index/$tool-index/$index.splicesites.txt
	chmod -R 777 /myvol1/index/${tool}-index/
	echo "Index is now saved under /$wd/index/$tool-index/$index"
}


### START here ############################################################################

# test filepaths
test_fasta
test_gtf

# Build Genome index if not already available
if ! test -f /$wd/index/$tool-index/$index/1.ht2; then build_index; fi

#make list of fastq files
mk_fastqlist

#make output directories
mk_outdir

### Start mapping ###

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
	#First attempt: Paired end mapping
	#...tag outputs with this flag to name it per fastqfile         "${line##*/}"

	hisat2 -q \
	-x /$wd/index/$tool-index/$index \
	-1 "$line"1.fastq -2 "$line"2.fastq \
	-S /$wd/$out/${line##*/}${tool}.sam \
	--known-splicesite-infile /$wd/index/$tool-index/$index.splicesites.txt

	#If paired end mapping fails, run unpaired mapping.
	trap 'second_attempt $line' ERR
done </$wd/tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
