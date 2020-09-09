#!/bin/bash

### Tool-specific variables ###
tool=contextmap

# use confic and function file
source /myvol1/config/mapping_config.sh
source /myvol1/func/mapping_func.sh

#To call contextmap as a command.
PATH="$PATH":/home/biodocker/ContextMap_Source_v2.7.9/:/home/biodocker/jre1.8.0_241/


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	#-reads 		A comma separated list of file paths to reads in fasta/fastq/zip/gz format. A single path for single-end mapping and two paths (#1 mates and #2 mates) for paired-end mapping
	#-aligner_name 	The following alignment tools are supported: "bwa", "bowtie1" or "bowtie2".
	#-aligner_bin 	The absolute path to the executable of the chosen aligner tool.
	#-indexer_bin 	The absolute path to the executable of the aligner's indexing tool (not needed for BWA).
	#-indices 		A comma separated list of paths to basenames of indices, which can be used by the chosen aligner.
	#-genome 		The path to a folder with reference sequences in fasta format (for each chromosome a separate file). It is important that the chromosome names contained in the index of the chosen unspliced alignment program (see 2. Obtaining ContextMap) are equal to the filenames of the files contained in the genome folder (for instance if the index contains a chromosome called "chr1", the genome folder must contain a file called "chr1.fa").
	#-o 			The path to the output directory.
	java -jar /home/biodocker/ContextMap_Source_v2.7.9/ContextMap_v2.7.9.jar \
		mapper \
		-reads ${line}?.fastq \
		-aligner_name bowtie2 \
		-o $out/${line##*/}${tool}.sam \
		-aligner_bin /home/biodocker/bin/\
		-indexer_bin /home/biodocker/bin/bowtie2-build \
		-indices $indexdir/$index \
		-genome $(ls $inputdir/$contextmap_fastadir)
}

#Build Genome index
build_index(){
	mkdir -p $indexdir/$index
	echo "compute index ..."
	for line in $fasta; do bowtie2-build -f $inputdir/$contextmap_fastadir/$line $indexdir/$line --threads $nthreads; done
	chmod -R 777 $indexdir
	echo "Index is now saved under $indexdir/$index"
}


### START here ############################################################################

# test filepaths
test_fasta

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$index.1.bt2; then build_index; fi fi

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

	#Parameters
	#-reads 		A comma separated list of file paths to reads in fasta/fastq/zip/gz format. A single path for single-end mapping and two paths (#1 mates and #2 mates) for paired-end mapping
	#-aligner_name 	The following alignment tools are supported: "bwa", "bowtie1" or "bowtie2".
	#-aligner_bin 	The absolute path to the executable of the chosen aligner tool.
	#-indexer_bin 	The absolute path to the executable of the aligner's indexing tool (not needed for BWA).
	#-indices 		A comma separated list of paths to basenames of indices, which can be used by the chosen aligner.
	#-genome 		The path to a folder with reference sequences in fasta format (for each chromosome a separate file). It is important that the chromosome names contained in the index of the chosen unspliced alignment program (see 2. Obtaining ContextMap) are equal to the filenames of the files contained in the genome folder (for instance if the index contains a chromosome called "chr1", the genome folder must contain a file called "chr1.fa").
	#-o 			The path to the output directory.
	
	java -jar /home/biodocker/ContextMap_Source_v2.7.9/ContextMap_v2.7.9.jar \
		mapper \
		-reads ${line}1.fastq,${line}2.fastq \
		-aligner_name bowtie2 \
		-o $out/${line##*/}${tool}.sam \
		-aligner_bin /home/biodocker/bin/bowtie2\
		-indexer_bin /home/biodocker/bin/bowtie2-build \
		-indices $indexdir/$index \
		-genome $(ls $inputdir/$contextmap_fastadir)

	#If paired end mapping fails, run unpaired mapping.
	trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner