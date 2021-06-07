#!/bin/bash

tool=kissplice
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

#cleaning up
trap cleaner EXIT


#make output directory
mk_outdir
#handle sam files
handlesamfiles $differential

#list all fastq files in fastqfolder and make list, seperating the file with " -r " to work with kissplice arguments; the first -r will be added manually
#if differential analysis: combine lists of casefastq-folder and controlfastq-folder into the fastqlist
#then run kissplice with this list of fastqfiles
if [ $differential = 0 ]
then
	fastqlist=$(ls -1p $fastqdir/*.fastq | xargs echo | sed 's/ / -r /g')
	echo Finding AS events...
	kissplice -r $fastqlist -o $outdir/kissplice
fi
if [ $differential = 1 ]
then
	casefastqlist=$(ls -1p $casefastq/*.fastq | xargs echo | sed 's/ / -r /g')
	controlfastqlist=$(ls -1p $controlfastq/*.fastq | xargs echo | sed 's/ / -r /g')
	echo Finding AS events...
	kissplice -r $casefastqlist -r $controlfastqlist -o $outdir/kissplice
fi


#aligning AS events back to reference genome
#check for STAR index
#if ! test -f $outdir/STAR_index/SAindex; then build_index; else echo STAR-index found, moving on...; fi
check_star_index

echo mapping AS event file back onto reference...
STAR --genomeDir $star_index --readFilesIn $outdir/kissplice/*_type_1.fa --outFileNamePrefix $outdir/ --runThreadN $ncores

echo classify events of KisSplice aligned to reference with kissplice2refgenome...
kissplice2refgenome -a $gtf -o $outdir/events --readLength $read_length $outdir/Aligned.out.sam


if [ $differential = 1 ]
then
	echo Starting differential analysis with kissDE...
	Rscript /docker_main/kissDE.R --counts $outdir/events.tsv --cores $ncores --caseprefix $caseprefix --controlprefix $controlprefix --casefastq $casefastq --controlfastq $controlfastq --out $outdir

fi
