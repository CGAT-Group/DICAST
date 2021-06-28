#!/bin/bash

### Tool-specific variables ###
tool=segemehl

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh

### logging ###
start_logging

### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	#Parameters
	# -d 	--database [<file>] 	list of filename(s) of fasta database sequence(s)
	# -q 	--query 	filename of query sequences
	# -p 	--mate 	filename of paired end sequences
	# -i 	--index 	filename of database index
	#X -j 	--index2 	filename of second database index
	#X -x 	--generate 	generate db index and store to disk
	#X -y 	--generate2 	generate second db index and store to disk
	#X -G 	--readgroupfile 	filename to read @RG header (default:none)
	#X -g 	--readgroupid 	read group id (default:none)
	# -t 	--threads 	start threads (default:1)
	#X -f 	--fullname 	write full fastq/a name to SAM
	segemehl.x \
		-d $fasta \
		-q ${line}?.fastq \
		--splits \
		-i $indexdir/$indexname \
		-t $ncores \
		-o $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}.sam
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index"
	segemehl.x -x $indexdir/$indexname -d $fasta
	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$indexname"
}

### START here ############################################################################

#cleaning up
trap cleaner EXIT

# test filepaths
test_fasta

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname; then build_index; fi fi

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

	# Parameters (X are not used)
	# -d 	--database [<file>] 	list of filename(s) of fasta database sequence(s)
	# -q 	--query 	filename of query sequences
	# -p 	--mate 	filename of paired end sequences
	# -i 	--index 	filename of database index
	#X -j 	--index2 	filename of second database index
	#X -x 	--generate 	generate db index and store to disk
	#X -y 	--generate2 	generate second db index and store to disk
	#X -G 	--readgroupfile 	filename to read @RG header (default:none)
	#X -g 	--readgroupid 	read group id (default:none)
	# -t 	--threads 	start threads (default:1)
	#X -f 	--fullname 	write full fastq/a name to SAM

	segemehl.x \
		-d $fasta \
		-q ${line}1.fastq  \
		-p ${line}2.fastq \
		--splits \
		-i $indexdir/$indexname \
		-t $ncores \
		-o $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}.sam

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

