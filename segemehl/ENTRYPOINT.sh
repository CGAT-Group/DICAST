#!/bin/bash

### Tool-specific variables ###
tool=segemehl

# use confic and function file
source /myvo11/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
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
		-d /$wd/$inputdir/$fasta \
		-q ${line}?.fastq \
		-i /$wd/index/$tool-index/$index \
		-t $nthreads \
		-o /$wd/$out/${line##*/}${tool}.sam
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p /$wd/index/$tool-index
	echo "compute index"
	segemehl.x -x /$wd/index/$tool-index/$index -d $(ls /$wd/$inputdir/$fasta)
	chmod -R 777 /$wd/index/$tool-index
	echo "Index is now saved under /$wd/index/$tool-index/$index.idx"
}

### START here ############################################################################

# test filepaths
test_fasta

# Build Genome index if not already available
if ! test -f /$wd/index/$tool-index/$index; then build_index; fi

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
	#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")
	
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
	-d /$wd/$inputdir/$fasta \
	-q ${line}1.fastq  \
	-p ${line}2.fastq \
	-i /myvol1/index/${tool}-index/$index \
	-t $nthreads \
	-o /$wd/$out/${line##*/}${tool}.sam

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </$wd/tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
