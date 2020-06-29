#!/bin/bash

# use config and function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=sgseq


### START here ###################

#test filepath for gtf and BAM+BAM-index
test_gtf $wd/$gtf
#test_bam $wd/$bamfolder/$bamfile



#make output directories and bamlist
mk_outdir $tool
readbamfiles

### Start AS event detection ###

echo compute ${tool} AS event detection...

for filename in $(cat $wd/$output/${tool:-unspecific}-output/bamlist)
do
	echo Starting SGSeq for $filename ...
        makebamfromsam $filename
	sample_out=$(mk_sample_out $filename)
	Rscript $wd/Rscripts/SGSeq.R --gtf $wd/$gtf --path_to_bam $filename --out $sample_out --cores $ncores
	wait
done


cleaner



