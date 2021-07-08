#!/bin/bash
tool=star

# use confic and function file
source /MOUNT/scripts/config.sh

source /MOUNT/scripts/mapping_func.sh

### logging ###
start_logging

#Build Genome index
build_index() {
	mkdir -p $star_index/$indexname
	echo "compute index ..."
	/docker_main/STAR-2.7.5c/bin/Linux_x86_64/STAR \
		--runMode genomeGenerate \
		--genomeDir $star_index/$indexname \
		--genomeFastaFiles $fasta \
		--runThreadN $ncores \
		--sjdbGTFfile $gtf \
		-sjdbOverhang 100 \
		--outFileNamePrefix $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}_${tool}
	chmod -R 777 $star_index
	echo "Index is now saved under $star_index/$indexname"
}

#Unpaired mapping command (EXPERIMENTAL)
second_attempt() {
for line1 in $(ls ${line}*.fastq| sed s/.fastq// );
do
	/docker_main/STAR-2.7.5c/bin/Linux_x86_64/STAR \
	--genomeDir $star_index/$indexname \
	--outFileNamePrefix $outdir/$(basename $(dirname $(dirname $line)))/${line1##*/}${tool} \
	--sjdbGTFfile $gtf  \
	--twopassMode Basic \
	--runThreadN $ncores \
	--outSAMstrandField intronMotif \
	--outSAMattributes NH HI AS nM NM XS \
	--readFilesIn ${line1}.fastq
done
}

### START here ############################################################################
#cleaning up
trap cleaner EXIT

# test filepaths
test_fasta
test_gtf

echo $gtf
#make output directories
#mk_outdir


# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $star_index/$indexname/genomeParameters.txt; then build_index; fi fi

#make list of fastq files
mk_fastqlist


### Start mapping ###

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
	trap 'second_attempt $line' ERR

	#First attempt: Paired end mapping
	echo mapping paired
	/docker_main/STAR-2.7.5c/bin/Linux_x86_64/STAR \
		--genomeDir $star_index/$indexname \
		--outFileNamePrefix $outdir/$(basename $(dirname $(dirname $line)))/${line##*/}${tool} \
		--sjdbGTFfile $gtf  \
		--twopassMode Basic \
		--runThreadN $ncores \
		--outSAMstrandField intronMotif \
		--outSAMattributes NH HI AS nM NM XS \
		--readFilesIn ${line}1.fastq ${line}2.fastq

	#If paired end mapping fails, run unpaired mapping. (EXPERIMENTAL)
done < /tmp/$tool-fastqlist

#wait for all processes to end
wait
