#!/bin/bash

#use config & function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=eventpointer

##### START here ###############


#test filepath for gtf and BAM+BAM-index
test_gtf $wd/$gtf
#test_bam $wd/$bamfolder/$bamfile


#make output directory
mk_outdir $tool

#handle SAM files
readsamfiles
for filename in $(cat $wd/$output/${tool:-unspecific}-output/samlist)
do
	makebamfromsam $filename
done



#### AS event detection mode ####
if [ $differential = 0 ]; then
	#make bamlist
	readbamfiles
	for filename in $(cat $wd/$output/${tool:-unspecific}-output/bamlist)
	do
	       	sample_out=$(mk_sample_out $filename)
	        echo Starting EventPointer for $filename in AS event detection mode...
        	Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores --out $sample_out --workdir $wd --bamfolder $bamfolder --bamfile $filename --differential $differential
        	wait
	done
cleaner
fi


#### DS mode ####
if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
	combined=$(combine_case_control $casefolder $controlfolder) #combine case and control folder into single folder
        Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores  --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $casefolder --controlfolder $controlfolder --differential $differential --combined $combined
	wait
	cleaner_diff
fi

