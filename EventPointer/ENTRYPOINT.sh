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
	echo Starting EventPointer in AS event detection mode...
	Rscript /docker_main/EventPointer.R --gtf $wd/$gtf --cores $ncores --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $wd/$bamfolder --differential $differential
	wait
	cleaner
fi


#### DS mode ####
if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
	combined=$(combine_case_control $casefolder $controlfolder) #combine case and control folder into single folder
        Rscript /docker_main/EventPointer.R --gtf $wd/$gtf --cores $ncores  --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $wd/$casefolder --controlfolder $wd/$controlfolder --differential $differential --combined $combined
	wait
	cleaner_diff
fi

