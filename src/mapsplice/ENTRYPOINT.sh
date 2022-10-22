#!/bin/bash

### Tool-specific variables ###
tool=mapsplice

# use confic and function file
source /MOUNT/scripts/config.sh

source /MOUNT/scripts/mapping_func.sh

### logging ###
start_logging

### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes (EXPERIMENTAL)
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	# Parameters
	# -a 		Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
	# -o FILE 	Output alignments to FILE [stdout].
	#-t INT 	Number of threads [3]. Minimap2 uses at most three threads when indexing target sequences,
	#			and uses up to INT+1 threads when mapping (the extra thread is for I/O, which is frequently
	#			idle and takes little CPU time).
	python /opt/conda/bin/mapsplice.py \
		-c $bowtie_fastadir \
		-x $indexdir/$indexname \
		--gene-gtf $gtf \
		-o $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool} \
		-p $ncores \
		-1 ${line}?.fastq
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "Build bowtie-index..."
	bowtie_fastadir_index=$(find $bowtie_fastadir -name "*" | sed -n 'H;${x;s/\n/,/g;s/^,//;p;}')

	# Usage: bowtie-build [options]* <reference_in> <ebwt_outfile_base>
	#    reference_in            comma-separated list of files with ref sequences
	#    ebwt_outfile_base       write Ebwt data to files with this dir/basename
	# Parameters:
	# --seed <int>            seed for random number generator

	bowtie-build \
		--seed 42 \
		#--threads $ncores \
		$bowtie_fastadir_index \
		$indexdir/$indexname
	chmod -R 777 $indexdir
}


### START here ############################################################################

#cleaning up
trap cleaner EXIT

# test filepaths
test_fasta
test_gtf

#update PATH
PATH=$PATH:/opt/conda/bin/mapsplice.py

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$indexname.rev.2.ebwt; then build_index; fi fi

#make list of fastq files
mk_fastqlist

#make output directories
#mk_outdir

### Start mapping ###

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
    mkdir -p $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool}

	#First attempt: Paired end mapping
	#...tag outputs with this flag to name it per fastqfile         "${line##*/}"

	# python mapsplice.py [options]* -c <Reference_Sequence> -x <Bowtie_Index> -1 <Read_List1> -2 <Read_List2>
	# Parameters (X means that the parameter is currently not used)
	# -c <string> 	The directory containing the sequence files of reference genome.
	# -x <string> 	The basename (including directory path) of Bowtie 1 index to be searched. The basename is the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / etc.
	# -p / --threads <int> 	Number of threads to be used for parallel aligning. Default is 1.
	# -o / --output <string> 	The directory in which MapSplice will write its output. Default is "./mapsplice_out/".
	# X --gene-gtf <string> 	Gene annotation file in GTF format, used to annotate fusion junctions

	python /opt/conda/bin/mapsplice.py \
		-c $bowtie_fastadir \
		-x $indexdir/$indexname \
		--gene-gtf $gtf \
		-o $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool} \
		-p $ncores \
		-1 ${line}1.fastq  \
		-2 ${line}2.fastq

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
