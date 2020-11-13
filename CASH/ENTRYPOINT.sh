#!/bin/bash

# use config and function file
tool=cash
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh


#test input-files
test_gtf $gtf

#making outdir
mk_outdir

#handle sam-files
handlesamfiles 1

#make bamlist for case & control
caselist=$(ls -1p $casebam/*.bam | xargs echo | sed 's/ /, /g')
controllist=$(ls -1p $controlbam/*.bam | xargs echo | sed 's/ /, /g')


echo Starting CASH run...

java -jar -Xmx10g /opt/cash_v2.2.1/cash.jar --Case:CASE $caselist --Control:CONTROL $controllist --GTF $gtf --Output $outdir

cleaner






