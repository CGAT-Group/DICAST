#!/bin/bash

# use config and function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=aspli


test_gtf $wd/$gtf
#test_bam $wd/$bamfolder/$bamfile

#making outdir

mk_outdir $tool

<<<<<<< HEAD
#checking for SAM files and building BAM-files with index with them
readsamfiles
for filename in $(cat $wd/$output/${tool:-unspecific}-output/samlist)
=======
#checking BAM files
readbamfiles
for filename in $(cat $wd/$output/${tool:-unspecific}-output/bamlist)
>>>>>>> 27968621b951552e19b0b60f0f2191487861456b
do
        makebamfromsam $filename
done

echo starting ASpli...

if [ $differential = 0 ]; then
	echo Starting ASpli in AS event detection mode...
	Rscript $wd/Rscripts/ASpli.R --gtf $wd/$gtf --cores $ncores --readLength $read_length --out $wd/$output/${tool:-unspecific}-output --workdir $wd --bamfolder $bamfolder --differential $differential
fi
if [ $differential = 1 ]; then
	echo Starting Aspli in DS analysis mode...
	Rscript $wd/Rscripts/ASpli.R --gtf $wd/$gtf --cores $ncores --readLength $read_length --out $wd/$output/${tool:-unspecific}-output --workdir $wd --casefolder $casefolder --controlfolder $controlfolder --differential $differential
fi


wait

cleaner 
