#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=spladder


#make output directory
mk_outdir $tool
out=$wd/$output/${tool:-unspecific}-output #local variable for output folder


#handle SAM file
readsamfiles
for filename in $(cat $out/samlist)
do
	makebamfromsam $filename
done


if [ $differential = 0 ]
then
	#list all .bam files in bamfolder with comma seperated
	bamlist=$(ls -1p $wd/$bamfolder/*.bam | xargs echo | sed 's/ /,/g')
	echo Starting spladder in event detection mode ...
	spladder build -b $bamlist -o $out -a $wd/$gtf --parallel $ncores -n $read_length
	cleaner
fi

if [ $differential = 1 ]
then
	caselist=$(ls -1p $wd/$casefolder/*.bam | xargs echo | sed 's/ /,/g')
	controllist=$(ls -1p $wd/$controlfolder/*.bam | xargs echo | sed 's/ /,/g')
	echo Starting spladder in DS mode ...
	echo "$caselist,$controllist"
	spladder build -b "$caselist,$controllist" -o $out -a $wd/$gtf --parallel $ncores -n $read_length

	echo testing for differential splicing  ...
	spladder test --conditionA $caselist --conditionB $controllist -o $out --parallel $ncores -n $read_length --labelA "CASE" --labelB "CONTROL"
	cleaner
fi
