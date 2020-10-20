#!/bin/bash

# use config and function file
tool=aspli
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh


#test input-files
test_gtf $gtf

#making outdir
mk_outdir $tool

#checking for SAM files and building BAM-files with index with them
readsamfiles
for filename in $(cat $outdir/samlist)
do
        makebamfromsam $filename
done

echo starting ASpli...

if [ $differential = 0 ]; then
	echo Starting ASpli in AS event detection mode...
<<<<<<< HEAD
	Rscript /docker_main/ASpli.R --gtf $wd/$gtf --cores $ncores --readLength $read_length --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $bamfolder --differential $differential
||||||| merged common ancestors
	Rscript $wd/Rscripts/ASpli.R --gtf $wd/$gtf --cores $ncores --readLength $read_length --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $bamfolder --differential $differential
=======
	Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir --bamfolder $bamdir --differential $differential
>>>>>>> unify-config
fi
if [ $differential = 1 ]; then
	echo Starting Aspli in DS analysis mode...
<<<<<<< HEAD
	Rscript /docker_main/ASpli.R --gtf $wd/$gtf --cores $ncores --readLength $read_length --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $casefolder --controlfolder $controlfolder --differential $differential
||||||| merged common ancestors
	Rscript $wd/Rscripts/ASpli.R --gtf $wd/$gtf --cores $ncores --readLength $read_length --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $casefolder --controlfolder $controlfolder --differential $differential
=======
	Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir --casefolder $casebam --controlfolder $controlbam --differential $differential
>>>>>>> unify-config
fi


wait

cleaner 
