#!/bin/bash

tool=spladder
source /MOUNT/scripts/config.sh 
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh


#make output directory
#mk_outdir $tool


#handle SAM file
handlesamfiles $differential

if [ $differential = 0 ]
then
	#list all .bam files in bamfolder with comma seperated
	bamlist=$(ls -1p $bamdir/*.bam | xargs echo | sed 's/ /,/g')
	echo Starting spladder in event detection mode ...
	spladder build -b $bamlist -o $outdir -a $gtf --parallel $ncores -n $read_length --output-txt-conf
	cleaner
fi

if [ $differential = 1 ]
then
	caselist=$(ls -1p $casebam/*.bam | xargs echo | sed 's/ /,/g')
	controllist=$(ls -1p $controlbam/*.bam | xargs echo | sed 's/ /,/g')
	echo Starting spladder in DS mode ...
	#echo "$caselist,$controllist"
	spladder build -b "$caselist,$controllist" -o $outdir -a $gtf --parallel $ncores -n $read_length --output-txt-conf

	echo testing for differential splicing  ...
	spladder test --conditionA $caselist --conditionB $controllist -o $outdir --parallel $ncores -n $read_length --labelA "CASE" --labelB "CONTROL"
	cleaner
fi
