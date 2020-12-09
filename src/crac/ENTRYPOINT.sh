#!/bin/bash

### Tool-specific variables ###
tool=crac

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	# Parameters
	# -i = index-folder/filename
	# -k = length of k-mers (k=22 recomended for human genome)
	# -r = name of fasta/fastq file containing reads
	# -o = output file
	# --detailed-sam = return a detailed sam file
	# -nb-threads = number of threads used
	# --stranded = strand specific rnaseq protocol
	crac \
		-i $indexdir/$indexname \
		-k 22 \
		-r ${line}?.fastq \
		-o $outdir/$(basename $(dirname $line))/${line##*/}${tool}.sam \
		--detailed-sam \
		--nb-threads $ncores
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index ..."
	crac-index index $indexdir/$indexname $(ls $fasta)
	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$index"
}


### START here ############################################################################

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname.ssa; then build_index; fi fi

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

	# Parameters
	# -i = index-folder/filename
	# -k = length of k-mers (k=22 recomended for human genome)
	# -r = name of fasta/fastq file containing reads
	# -o = output file
	# --detailed-sam = return a detailed sam file
	# -nb-threads = number of threads used
	# --stranded = strand specific rnaseq protocol

	crac \
		-i $indexdir/$indexname \
		-k 22 \
		-r ${line}1.fastq ${line}2.fastq \
		-o $outdir/$(basename $(dirname $line))/${line##*/}${tool}.sam \
		--detailed-sam \
		--nb-threads $ncores \
		--stranded

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done < /tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
