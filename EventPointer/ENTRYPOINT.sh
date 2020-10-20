#!/bin/bash

tool=eventpointer
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh


##### START here ###############


#test filepath for gtf and BAM+BAM-index
test_gtf $gtf

#make output directory
mk_outdir $tool

#handle SAM files
readsamfiles
for filename in $(cat $outdir/samlist)
do
	makebamfromsam $filename
done



#### AS event detection mode ####
if [ $differential = 0 ]; then
	echo Starting EventPointer in AS event detection mode...
<<<<<<< HEAD
	Rscript /docker_main/EventPointer.R --gtf $wd/$gtf --cores $ncores --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $wd/$bamfolder --differential $differential
||||||| merged common ancestors
	Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $wd/$bamfolder --differential $differential
=======
	Rscript /docker_main/EventPointer.R --gtf $gtf --cores $ncores --out $outdir --bamfolder $bamdir --differential $differential
>>>>>>> unify-config
	wait
	cleaner
fi


#### DS mode ####
if [ $differential = 1 ]; then
        echo Starting EventPointer in DS analysis mode...
<<<<<<< HEAD
	combined=$(combine_case_control $casefolder $controlfolder) #combine case and control folder into single folder
        Rscript /docker_main/EventPointer.R --gtf $wd/$gtf --cores $ncores  --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $wd/$casefolder --controlfolder $wd/$controlfolder --differential $differential --combined $combined
||||||| merged common ancestors
	combined=$(combine_case_control $casefolder $controlfolder) #combine case and control folder into single folder
        Rscript $wd/Rscripts/EventPointer.R --gtf $wd/$gtf --cores $ncores  --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $wd/$casefolder --controlfolder $wd/$controlfolder --differential $differential --combined $combined
=======
	combined=$(combine_case_control $casebam $controlbam) #combine case and control folder into single folder
        Rscript /docker_main/EventPointer.R --gtf $gtf --cores $ncores  --out $outdir --casefolder $casefolder --controlfolder $controlfolder --differential $differential --combined $combined
>>>>>>> unify-config
	wait
	cleaner_diff
fi

