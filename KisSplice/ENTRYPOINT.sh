#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=kissplice


build_index(){
	echo "Did not find STAR index in $out/STAR_index/ ; building it now..."
	STAR --runMode genomeGenerate --genomeDir $out/tmp/STAR_index --genomeFastaFiles $wd/$fasta --sjdbGTFfile $wd/$gtf --runThreadN $ncores --outFileNamePrefix $out/tmp/Star_mapped_ --sjdbOverhang 100
	echo "STAR index built and saved to $out/STAR_index/"
}

#make output directory
mk_outdir $tool
out=$wd/$output/${tool:-unspecific}-output #local variable for output folder


#list all fastq files in fastqfolder and make list, seperating the file with " -r " to work with kissplice arguments; the first -r will be added manually
#if differential analysis: combine lists of casefastq-folder and controlfastq-folder into the fastqlist
#then run kissplice with this list of fastqfiles
if [ $differential = 0 ]
then
	fastqlist=$(ls -1p $wd/$fastqfolder/*.fastq | xargs echo | sed 's/ / -r /g')
	echo Finding AS events...
	kissplice -r $fastqlist -o $out/tmp
fi 
if [ $differential = 1 ]
then
	casefastqlist=$(ls -1p $wd/$casefastq/*.fastq | xargs echo | sed 's/ / -r /g')
	controlfastqlist=$(ls -1p $wd/$controlfastq/*.fastq | xargs echo | sed 's/ / -r /g')
	echo Finding AS events...
	kissplice -r $casefastqlist -r $controlfastqlist -o $out/tmp
fi	


#aligning AS events back to reference genome
#check for STAR index
if ! test -f $out/tmp/STAR_index/SAindex; then build_index; else echo STAR-index found, moving on...; fi

echo mapping AS event file back onto reference...
STAR --genomeDir $out/tmp/STAR_index --readFilesIn $out/tmp/*_type_1.fa --outFileNamePrefix $out/ --runThreadN $ncores

echo classify events of KisSplice aligned to reference with kissplice2refgenome...
kissplice2refgenome -a $wd/$gtf -o $out/events --readLength $read_length $out/Aligned.out.sam 


if [ $differential = 1 ]
then
	echo Starting differential analysis with kissDE...
	Rscript $wd/Rscripts/kissDE.R --counts $out/events.tsv --cores $ncores --caseprefix $caseprefix --controlprefix $controlprefix --casefastq $wd/$casefastq --controlfastq $wd/$controlfastq --out $out

fi


#cleaner


