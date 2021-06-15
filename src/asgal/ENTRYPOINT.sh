#!/bin/bash

tool=asgal
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

#cleaning up
trap cleaner EXIT

#make output directory
mk_outdir
handlesamfiles
#one asgal run needs a pair of fastq-files; in the config file the user gave the suffixes which correspond to the partnered fastq-files
#save all partners with same suffix in array
partner1fastqlist=($(ls -1p $fastqdir/*$fastqpair1suffix | xargs echo | sed 's/ / /g' | uniq))
partner2fastqlist=($(ls -1p $fastqdir/*$fastqpair2suffix | xargs echo | sed 's/ / /g' | uniq))

#check if both arrays are of same size
if [ ${#partner1fastqlist[@]} != ${#partner2fastqlist[@]} ]
then
	echo Did not find a pair for each fastq-file! Please check that each fastq file has a forward and reverse version and they share the suffixed, given in the config file.
	exit 1
fi

#get number of partners 
nPartners=${#partner2fastqlist[@]}

#iterate over fastqlists and run asgal for each fastq-pair
for ((i=0;i<nPartners;++i)); do
	#get the two corresponding fastq files from array
	fastq1=${partner1fastqlist[i]}
	fastq2=${partner2fastqlist[i]}
	echo "Starting ASGAL run for $fastq1 and $fastq2 ..."
	#create output folder for fastq-pair (named by first file)
	sample_out=$(mk_sample_out $fastq1)
	#run ASGAL
	/docker_main/galig/asgal --multi -g $fasta -a $gtf -t $transcript -s $fastq1 -s2 $fastq2 -o $sample_out -@ $ncores
	wait
done

