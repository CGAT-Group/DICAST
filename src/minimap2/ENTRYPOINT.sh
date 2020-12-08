#!/bin/bash

### Tool-specific variables ###
tool=minimap2

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
	# -a 		Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default. 
	# -o FILE 	Output alignments to FILE  [stdout]. 
	# -t INT 	Number of threads [3]. Minimap2 uses at most three threads when indexing target sequences, 
	#			and uses up to INT+1 threads when mapping (the extra thread is for I/O, which is frequently 
	#			idle and takes little CPU time). 
	minimap2 \
		-a \
		-t $ncores \
		-o $outdir/$(basename $(dirname $line))/${line##*/}$tool.sam \
		$indexdir/$indexname \
		$line?.fastq
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index ..."
	minimap2 -d $indexdir/$indexname $fasta
	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$index"
}

### START here ############################################################################

# test filepaths
test_fasta
	
# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname; then build_index; fi fi

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
	# -a 		Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default. 
	# -o FILE 	Output alignments to FILE [stdout]. 
	#-t INT 	Number of threads [3]. Minimap2 uses at most three threads when indexing target sequences, 
	#			and uses up to INT+1 threads when mapping (the extra thread is for I/O, which is frequently 
	#			idle and takes little CPU time). 

	minimap2 \
		-a \
		-t $ncores \
		-o $outdir/$(basename $(dirname $line))/${line##*/}${tool}.sam \
		$indexdir/$indexname \
		${line}1.fastq \
		${line}2.fastq

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
