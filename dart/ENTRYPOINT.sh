#!/bin/bash

### Tool-specific variables ###
tool=dart

# use confic and function file
source /myvo11/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	dart \
		-i /$wd/index/$tool-index/$index \
		-f $line?.fastq  \
		-o /$wd/$out/${line##*/}$tool.sam
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p /$wd/index/$tool-index
	echo "compute index ..."
	bwt_index $(ls /$wd/$fasta) /$wd/index/$tool-index/$index
	chmod -R 777 /$wd/index/$tool-index
	echo "Index is now saved at /$wd/index/$tool-index/$index"
}

### START here ############################################################################

# test filepaths
test_fasta

# Build Genome index if not already available
if ! test -f /$wd/index/$tool-index/$index.sa; then build_index; fi

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
	#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")

	dart \
		-i /$wd/index/$tool-index/$index \
		-f ${line}1.fastq \
		-f2 ${line}2.fastq \
		-o /$wd/$out/${line##*/}${tool}.sam
		
	#If paired end mapping fails, run unpaired mapping.
	trap 'second_attempt $line' ERR
done </$wd/tmp/$tool-fastqlist


# wait for all processes to end
wait
cleaner
