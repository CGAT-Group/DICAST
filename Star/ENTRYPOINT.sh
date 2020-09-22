#!/bin/bash
tool=star
star_version=2.7.5c
source /myvol1/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


#Build Genome index
build_index() {
	mkdir $indexdir; chmod 777 $indexdir
	/docker_main/STAR-${star_version}/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir $indexdir --genomeFastaFiles $inputdir/$fasta --runThreadN $nthreads --sjdbGTFfile $inputdir/$gtf -sjdbOverhang 100 --outFileNamePrefix $out/Star_mapped_${line##*/}
}

#Unpaired mapping command
second_attempt() {
for line1 in $(ls ${line}*.fastq| sed s/.fastq// );
do
	/docker_main/STAR-${star_version}/bin/Linux_x86_64/STAR --genomeDir $indexdir --outFileNamePrefix $out/Star_mapped_${line1##*/} --sjdbGTFfile $inputdir/$gtf  --twopassMode Basic --runThreadN $nthreads --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn ${line1}.fastq
done
}


#START here: Make a list of fastq files
mk_fastqlist


#make output directories
mk_outdir

#test filepaths for fasta and indexing
test_fasta
test_gtf

mk_outdir
if ! test -f $indexdir/genomeParameters.txt; then build_index; fi


#Iterate list with paired end map command first
while read -r line; do
	trap 'second_attempt $line' ERR

	#First attempt: Paired end mapping
echo mapping paired
	/docker_main/STAR-${star_version}/bin/Linux_x86_64/STAR --genomeDir  $indexdir --outFileNamePrefix $out/Star_mapped_${line##*/} --sjdbGTFfile $inputdir/$gtf  --twopassMode Basic --runThreadN $nthreads --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn ${line}1.fastq ${line}2.fastq

	#If paired end mapping fails, run unpaired mapping.
done < /tmp/$tool-fastqlist

#wait for all processes to end
wait

#cleaning up
trap cleaner EXIT
