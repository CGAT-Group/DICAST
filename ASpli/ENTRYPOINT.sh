#!/bin/bash

# use config and function file
tool=aspli
source template/config.sh
source template/asevent_config.sh
source template/asevent_func.sh


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
	Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir --workdir $workdir --bamfolder $bamdir --differential $differential
fi
if [ $differential = 1 ]; then
	echo Starting Aspli in DS analysis mode...
	Rscript /docker_main/ASpli.R --gtf $gtf --cores $ncores --readLength $read_length --out $outdir --workdir $workdir --casefolder $casefolder --controlfolder $controlfolder --differential $differential
fi


wait

cleaner 
