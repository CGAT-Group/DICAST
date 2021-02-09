#!/bin/bash

# use config and function file
tool=edger
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh


#test input-files
test_gtf $gtf

#making outdir
#mk_outdir $tool

#checking for SAM files and building BAM-files with index with them
handlesamfiles 1


echo starting edgeR...

Rscript /docker_main/edgeR.R --gtf $gtf --casedir $casebam --controldir $controlbam --output $outdir
wait

cleaner



