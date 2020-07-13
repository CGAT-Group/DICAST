#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=asgal

#make output directory
mk_outdir $tool


#one asgal run needs a pair of fastq-files; in the config file the user gave the suffixes which correspond to the partnered fastq-files
#save all partners with same suffix in array
partner1fastqlist=($(ls -1p $wd/$fastqfolder/*$fastqpair1suffix | xargs echo | sed 's/ / /g' | uniq))
partner2fastqlist=($(ls -1p $wd/$fastqfolder/*$fastqpair2suffix | xargs echo | sed 's/ / /g' | uniq))

#check if both arrays are of same size
if [ ${#partner1fastqlist[@]} != ${#partner2fastqlist[@]} ]
then
	echo Did not find a pair for each fastq-file! Please check that each fastq file has a forward and reverse version and they share the suffixed, given in the config file.
	exit 1
fi

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
	/galig/asgal --multi -g $wd/$fasta -a $wd/$gtf -t $wd/$transcript -s $fastq1 -s2 $fastq2 -o $sample_out -@ $ncores
	wait
done

cleaner
