#!/bin/bash
tool=star

source /myvol1/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


#Build Genome index
build_index() {
	/docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir /$wd/$out/genomedir --genomeFastaFiles /$wd/$fasta --runThreadN $nthreads --sjdbGTFfile /$wd/$gtf -sjdbOverhang 100 --outFileNamePrefix /$wd/$out/Star_mapped_${line##*/}
}

#Unpaired mapping command
second_attempt() {
	/docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir /$wd/$out/genomedir --outFileNamePrefix /$wd/$out/Star_mapped_${line##*/} --sjdbGTFfile /$wd/$gtf  --twopassMode Basic --runThreadN $nthreads --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn ${line}?.fastq

}


#START here: Make a list of fastq files
mk_fastqlist


#make output directories
mk_outdir

#test filepaths for fasta and indexing
test_fasta
test_gtf
mkdir /$wd/$out/genomedir
if ! test -f /$wd/$out/genomedir/SAindex; then build_index; fi


#Iterate list with paired end map command first
while read -r line; do

	#First attempt: Paired end mapping
echo mapping paired
	/docker_main/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir  /$wd/$out/genomedir --outFileNamePrefix /$wd/$out/Star_mapped_${line##*/} --sjdbGTFfile /$wd/$gtf  --twopassMode Basic --runThreadN $nthreads --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM XS --readFilesIn ${line}1.fastq ${line}2.fastq

	#If paired end mapping fails, run unpaired mapping.
	trap 'second_attempt $line' ERR
done < /$wd/tmp/$tool-fastqlist

#wait for all processes to end
wait

#cleaning up
trap cleaner EXIT
