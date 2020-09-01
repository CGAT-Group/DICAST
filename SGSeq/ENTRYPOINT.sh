#!/bin/bash

# use config and function file
source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=sgseq


### START here ###################

#tests
mk_outdir $tool
test_gtf $wd/$gtf
out=$wd/$output/${tool:-unspecific}-output #local variable for output folder


#handle sam files
readsamfiles
for filename in $(cat $wd/$output/${tool:-unspecific}-output/samlist)
do
        makebamfromsam $filename
done

#build bamlist
readbamfiles

### Start AS event detection ###

echo compute ${tool} AS event detection...

for filename in $(cat $wd/$output/${tool:-unspecific}-output/bamlist)
do
	echo Starting SGSeq for $filename ...
	sample_out=$(mk_sample_out $filename)
	Rscript $wd/Rscripts/SGSeq.R --gtf $wd/$gtf --path_to_bam $filename --out $sample_out --cores $ncores
	wait
done


cleaner



