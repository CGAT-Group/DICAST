#!/bin/bash

### Tool-specific variables ###
tool=hisat

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh

### logging ###
start_logging

#Update PATH
PATH=$PATH:/docker_main/hisat2-2.2.1

### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	hisat2 -q \
	-x $indexdir/$indexname \
	-U "$line"?.fastq \
	-S $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}.sam \
	--known-splicesite-infile $indexdir/${gtfname}_splicesites.txt
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index ..."
	hisat2-build  $fasta $indexdir/$indexname
	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$indexname"
}

extract_splice_sites() {
	python3 /docker_main/hisat2-2.2.1/hisat2_extract_splice_sites.py $gtf > $indexdir/${gtfname}_splicesites.txt
	chmod -R 777 $indexdir
	echo "Extracted plicesites are now saved under $indexdir/${gtfname}_splicesites.txt"
}

### START here ############################################################################

#cleaning up
trap cleaner EXIT

# test filepaths
test_fasta
test_gtf

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/${indexname}.4.ht2; then build_index; fi fi

if $recompute_index; then extract_splice_sites; else if ! test -f $indexdir/${gtfname}_splicesites.txt; then extract_splice_sites; fi fi


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

	hisat2 -q \
	--known-splicesite-infile $indexdir/${gtfname}_splicesites.txt \
	-x $indexdir/$indexname \
	-1 "$line"1.fastq -2 "$line"2.fastq \
	-S $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}.sam


	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
