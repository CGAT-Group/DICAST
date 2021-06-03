#!/bin/bash

tool=irfinder
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

#cleaning up
trap cleaner EXIT


mk_outdir
#######################################
### build IRFinder custom reference ###
#######################################

test_gtf $gtf
test_fasta $fasta

#handle SAM files
handlesamfiles 0

#link/move gtf and fasta files into output folder
echo linking annotation files into reference folder...
mkdir -p $outdir/irfinder_index

#check if GTF is for IRFinder, else attempt to 'fix' issue
if grep biotype $gtf ; then
	ln -sf $gtf $outdir/irfinder_index/transcripts.gtf
else
	echo Attempting to \'fix\' gtf by adding \'biotype\' like IRFinder wants...
	python /docker_main/gtf_for_irfinder.py $gtf $outdir/irfinder_index/transcripts.gtf
fi

ln -sf $fasta $outdir/irfinder_index/genome.fa


#build reference
echo building reference...
IRFinder -m BuildRefProcess -r $outdir/irfinder_index 
wait
echo reference built, moving on...



#####################
#### quantify IR ####
#####################

#readfastqs $fastqfolder
#fastqs=$(cat $out/fastqlist)

echo looking for intron retention events...
if [[ $use_bam_input_files -eq 0 ]]
then
	#one irfinder run needs a pair of fastq-files; in the config file the user gave the suffixes which correspond to the partnered fastq-files
	#save all partners with same suffix in array
	partner1fastqlist=($(ls -1p $fastqdir/*$fastqpair1suffix | xargs echo | sed 's/ / /g' | uniq))
	partner2fastqlist=($(ls -1p $fastqdir/*$fastqpair2suffix | xargs echo | sed 's/ / /g' | uniq))


	#check if both arrays are of same size
	if [ ${#partner1fastqlist[@]} != ${#partner2fastqlist[@]} ] 
	then 	
		echo Did not find a pair for each fastq-file! Please check that each fastq file has a forward and reverse version and they share the suffixed, given in the config file. 	
		exit 1 
	fi

	
	nPartners=${#partner2fastqlist[@]}
	#iterate over fastqlists and run irfinder for each fastq-pair
	for ((i=0;i<nPartners;++i)); do
		#get the two corresponding fastq files from array
		fastq1=${partner1fastqlist[i]}
		fastq2=${partner2fastqlist[i]}
		echo "Starting IRFinder run for $fastq1 and $fastq2 ..."
		#create output folder for fastq-pair (named by first file)
		sample_out=$(mk_sample_out $fastq1)
		#run irdinder
		IRFinder -r $outdir/irfinder_index -d $sample_out $fastq1 $fastq2
		wait
	done



elif [[ $use_bam_input_files -eq 1  ]]
then
	readbamfiles
	bams=$(cat /tmp/controlbamlist)
	for bam in $bams
	do
       		IRFinder -m BAM -r $outdir/irfinder_index -d $outdir $bam
       		wait
	done
fi