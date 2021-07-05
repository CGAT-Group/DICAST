#!/bin/bash

tool=asgal
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

#cleaning up
trap cleaner EXIT


#make sure the $transcript file exists
if [[ ! -s $transcript ]]
then gffread $gff -g $fasta -w $transcript; new transcript file generated at $transcript
	else for i in $(seq 1 9999) 
		do if ! ls ${transcript}$i 2>/dev/null 
			then mv $transcript ${transcript}$i && break
		fi
	done ;	echo older transcript file found at $transcript, renamed to: ${transcript}$i
fi


#make output directory
#mk_outdir
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
	fastqname=$(basename $fastq1 .fastq | sed 's/..$//')
	sample_out=$(mk_sample_out $fastqname)
	unified_outdir_name="${sample_out}_${tool}_unified"

	#run ASGAL
	/galig/asgal --multi -g $fasta -a $gtf -t $transcript -s $fastq1 -s2 $fastq2 -o $sample_out -@ $ncores --allevents
	cat $sample_out/ASGAL/*.events.csv | head -n1 > $sample_out/ASGAL/all.events.csv
	cat $sample_out/ASGAL/*.events.csv | grep -v 'Type,Start' >> $sample_out/ASGAL/all.events.csv
	mkdir -p $unified_outdir_name
	if [ $combine_events = 0 ]; 
		then
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -a $sample_out/ASGAL/all.events.csv -out $unified_outdir_name -gtf $gtf
		else
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -a $sample_out/ASGAL/all.events.csv -out $unified_outdir_name -gtf $gtf -comb
		fi
done
