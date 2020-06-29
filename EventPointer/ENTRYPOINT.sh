#!/bin/bash

#use config & function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=eventpointer

##### START here ###############


#test filepath for gtf and BAM+BAM-index
test_gtf $wd/$gtf
#test_bam $wd/$bamfolder/$bamfile


#make output directory and check BAM + BAM index
mk_outdir $tool
readbamfiles


#### AS event detection mode ####
for filename in $(cat $wd/$output/${tool:-unspecific}-output/bamlist)
do
	if [ $differential = 0 ]; then
		makebamfromsam $filename
        	sample_out=$(mk_sample_out $filename)
	        echo Starting EventPointer for $filename in AS event detection mode...
        	Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores --out $sample_out --workdir $wd --bamfolder $bamfolder --bamfile $filename --differential $differential
        	wait
        	rm -f $wd/$output/${tool:-unspecific}-output/ASpli_binFeatures.log
	fi
done
cleaner


#### DS mode ####
if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
	combined=$(combine_case_control) #combine case and control folder into single folder
        Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores  --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $casefolder --controlfolder $controlfolder --differential $differential --combined $combined
	wait
	rm -f $wd/$output/${tool:-unspecific}-output/ASpli_binFeatures.log
	cleaner_diff
fi

