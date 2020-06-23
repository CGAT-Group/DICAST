#!/bin/bash

### Tool-specific variables ###
tool=subjunc

# use confic and function file
source /myvo11/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	# Parameters
	# -r	Give the name of an input file (multiple files are allowed tobe provided toalignandsubjuncfunctions inRsubread).
	#		For paired-end read data, this gives the first read file andthe other read file should be provided via the -R option.
	#		Supported input formats include FASTQ/FASTA (uncom-pressed or gzip compressed)(default), SAM and BAM.
	# -R	Provide name of the second read file from paired-end data.The program will switch to paired-end read mapping modeif  this  file  is  provided.   
	#		(multiple  files  are  allowed  to  beprovided toalignandsubjuncfunctions inRsubread)
	# -i	Specify  the  base  name  of  the  index.   The  index  used  by sublong aligner must be a full index and has only one block,ie.  
	#		‘-F’ and ‘-B’ options must be specified when buildingindex with subread-buildindex.
	# -T	Specify  the  number  of  threads/CPUs  used  for  mapping.The value should be between 1 and 32.  1 by default
	# --SAMoutput	Specify that mapping results are saved into a SAM formatfile.
	# -o	Give the name of output file.  The default output format isBAM. All reads are included in mapping output, includingboth  mapped  and  unmapped  reads,  
	#		and  they  are  in  the same order as in the input file.

	subjunc \
		-r ${line}?.fastq \
		-i /$wd/index/$tool-index/$index \
		-T $nthreads \
		−−SAMoutput \
		-o /$wd/$out/${line##*/}$tool.sam
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p /$wd/index/$tool-index
	echo "compute index ..."
	subread-buildindex $(ls /$wd/$inputdir/$fasta) -o /$wd/index/$tool-index/$index -F -B 
	chmod -R 777 /$wd/index/$tool-index
	echo "Index is now saved under /$wd/index/$tool-index/$index"
}

### START here ############################################################################

# test filepaths
test_fasta
	
#Update PATH variable
PATH=$PATH:/opt/subread-2.0.0-Linux-x86_64/bin/ 

# Build Genome index if not already available
if ! test -f /$wd/index/$tool-index/$index.reads; then build_index; fi

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

# -r	Give the name of an input file (multiple files are allowed tobe provided toalignandsubjuncfunctions inRsubread).
#		For paired-end read data, this gives the first read file andthe other read file should be provided via the -R option.
#		Supported input formats include FASTQ/FASTA (uncom-pressed or gzip compressed)(default), SAM and BAM.
# -R	Provide name of the second read file from paired-end data.The program will switch to paired-end read mapping modeif  this  file  is  provided.   
#		(multiple  files  are  allowed  to  beprovided toalignandsubjuncfunctions inRsubread)
# -i	Specify  the  base  name  of  the  index.   The  index  used  by sublong aligner must be a full index and has only one block,ie.  
#		‘-F’ and ‘-B’ options must be specified when buildingindex with subread-buildindex.
# -T	Specify  the  number  of  threads/CPUs  used  for  mapping.The value should be between 1 and 32.  1 by default
# --SAMoutput	Specify that mapping results are saved into a SAM formatfile.
# -o	Give the name of output file.  The default output format isBAM. All reads are included in mapping output, includingboth  mapped  and  unmapped  reads,  
#		and  they  are  in  the same order as in the input file.


subjunc \
	-i /$wd/index/$tool-index/$index \
	-r ${line}1.fastq \
	-R ${line}2.fastq \
	-o /$wd/$out/${line##*/}${tool}.sam \
	-T $nthreads \
	−−SAMoutput

	#If paired end mapping fails, run unpaired mapping.
	trap 'second_attempt $line' ERR
done </$wd/tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner

