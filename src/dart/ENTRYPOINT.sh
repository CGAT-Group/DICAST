#!/bin/bash

### Tool-specific variables ###
tool=dart

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	dart \
		-i $indexdir/$indexname \
		-f $line?.fastq  \
		-o $outdir/$(basename $(dirname $line))/${line##*/}$tool.sam
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index ..."
	bwt_index $fasta $indexdir/$indexname
	chmod -R 777 $indexdir
	echo "Index is now saved at $indexdir/$indexname"
}

### START here ############################################################################

# test filepaths
test_fasta

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname.sa; then build_index; fi fi

#make list of fastq files
mk_fastqlist

#make output directories
#mk_outdir

### Start mapping ###

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
	#First attempt: Paired end mapping
	#...tag outputs with this flag to name it per fastqfile         "${line##*/}"

	dart \
		-i $indexdir/$indexname \
		-f ${line}1.fastq \
		-f2 ${line}2.fastq \
		-o $outdir/$(basename $(dirname $line))/${line##*/}${tool}.sam
		
	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist


# wait for all processes to end
wait
cleaner
