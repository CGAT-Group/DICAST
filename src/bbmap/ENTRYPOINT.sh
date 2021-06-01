#!/bin/bash

### Tool-specific variables ###
tool=bbmap

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
 	/bbmap/bbmap.sh \
		in=${line}?.fastq \
		outm=$outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}.sam \
		outu=$outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}_unmapped.sam
}

build_index() {
	mkdir -p $indexdir/$indexname
	echo "compute index ..."
	/bbmap/bbmap.sh ref=$fasta path=$indexdir/$indexname overwrite=false
	chmod -R 777 $indexdir
	echo "Index is now saved at $indexdir/$indexname"
}

### START here ############################################################################

#cleaning up
trap cleaner EXIT

# test filepaths
test_fasta

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -d $indexdir/$indexname; then build_index; fi fi

#make list of fastq files
mk_fastqlist

#make output directories
#mk_outdir

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
	#First attempt: Paired end mapping
	#...tag outputs with this flag to name it per fastqfile         "${line##*/}"

	/bbmap/bbmap.sh \
		in=${line}1.fastq \
		in2=${line}2.fastq \
		ref=$fasta \
		path=$indexdir/$indexname \
		intronlen=20 \
		xstag=us \
		outm=$outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}.sam \
		outu=$outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}_unmapped.sam

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait