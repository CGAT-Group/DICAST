#!/bin/bash

### Tool-specific variables ###
tool=gsnap

# use confic and function file
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/mapping_config.sh
source /MOUNT/scripts/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	# Usage: gsnap [OPTIONS...] <FASTA file>, or
	# Parameters (X means that the parameter is currently not used)
	# --db	Genome database
	# --dir Genome directory
	# --nthreads	Number of worker threads
	# --output-file File name for a single stream of output results.
	# --format=STRING Another format type, other than default. Currently implemented: sam, m8 (BLAST tabular format)
	# X --ordered      Print output in same order as input (relevant only if there is more than one worker thread)
	# ...
	gsnap \
		--db $indexname \
		--dir $indexdir \
		--output-file $outdir/${line##*/}${tool}.sam \
		--format sam \
		--nthreads $ncores \
		${line}?.fastq

}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index ..."
	# Usage: gmap_build [options...] -d <genomename> <fasta_files>
	# -D, --dir=STRING          Destination directory for installation (defaults to gmapdb directory specified at configure time)
	# -d, --db=STRING           Genome name

	gmap_build \
	-db $indexname \
	-dir $indexdir \
	$fasta

	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$indexname"
}

### START here ############################################################################

# test filepaths
test_fasta

# Build genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname/$indexname.contig; then build_index; fi fi

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
 
	# Usage: gsnap [OPTIONS...] <FASTA file>
	# Parameters (X means that the parameter is currently not used)
	# --db	Genome database
	# --dir Genome directory
	# --nthreads	Number of worker threads
	# --output-file File name for a single stream of output results.
	# --format=STRING Another format type, other than default. Currently implemented: sam, m8 (BLAST tabular format)
	# X --ordered      Print output in same order as input (relevant only if there is more than one worker thread)
	# ...

	gsnap \
	--db $indexname \
	--dir $indexdir \
	--output-file $outdir/${line##*/}${tool}.sam \
	--format sam \
	--nthreads $ncores \
	${line}1.fastq ${line}2.fastq

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
