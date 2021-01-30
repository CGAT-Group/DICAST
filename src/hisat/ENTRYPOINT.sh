#!/bin/bash

### Tool-specific variables ###
tool=hisat

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh

#Update PATH
PATH=$PATH:/docker_main/hisat2-2.0.0-beta

### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	hisat2 -q \
	-x $indexdir/$indexname \
	-U "$line"?.fastq \
	-S $outdir/$(basename $(dirname $line))/${line##*/}${tool}.sam \
	--known-splicesite-infile $indexdir/$index.splicesites.txt
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index ..."
	hisat2-build  $fasta $indexdir/$indexname && python /docker_main/hisat2-2.0.0-beta/extract_splice_sites.py $gtf > $indexdir/$indexname/splicesites.txt
	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$indexname"
}


### START here ############################################################################

# test filepaths
test_fasta
test_gtf

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname_1.ht2; then build_index; fi fi

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
	-x $indexdir/$indexname \
	-1 "$line"1.fastq -2 "$line"2.fastq \
	-S $outdir/$(basename $(dirname $line))/${line##*/}${tool}.sam \
	--known-splicesite-infile $indexdir/$indexname/splicesites.txt

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
