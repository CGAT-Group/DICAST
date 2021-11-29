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
if [[ -s $transcript ]]
then
	for i in $(seq 1 9999) 
		do if ! ls ${transcript}$i 2>/dev/null 
			then mv $transcript ${transcript}$i && break
		fi
	done ;	echo older transcript file found at $transcript, renamed to: ${transcript}$i
fi
/galig/gffread/gffread $gff -g $fasta -w $transcript && echo new transcript file generated at: ls transcript$ $(ls -lh $transcript)


#make output directory
#mk_outdir
handlesamfiles 2>/dev/null
#one asgal run needs a pair of fastq-files; in the config file the user gave the suffixes which correspond to the partnered fastq-files
#save all partners with same suffix in array
partner1fastqlist=($(ls -1p $fastqdir/*$fastqpair1suffix | uniq))
partner2fastqlist=($(ls -1p $fastqdir/*$fastqpair2suffix | uniq))

#check if both arrays are of same size
if [ ${#partner1fastqlist[@]} != ${#partner2fastqlist[@]} ]
then
	echo Did not find a pair for each fastq-file! Please check that each fastq file has a forward and reverse version and they share the suffixed, given in the config file.
	exit 1
fi

#get number of partners
nPartners=${#partner2fastqlist[@]}
echo fastqlist: ${partner1fastqlist[*]}

#iterate over fastqlists and run asgal for each fastq-pair
for ((i=0;i<nPartners;++i)); do
	#get the two corresponding fastq files from array
	fastq1=${partner1fastqlist[i]}
	fastq2=${partner2fastqlist[i]}
	echo "Starting ASGAL run for $fastq1 and $fastq2 ..."
	#create output folder for fastq-pair (named by first file)
	fastqname=$(basename $fastq1 .fastq | sed 's/..$//')
	sample_out=$(mk_sample_out $fastqname)


	#run ASGAL
	/galig/asgal --multi -g $fasta -a $gtf -t $transcript -s $fastq1 -s2 $fastq2 -o $sample_out -@ $ncores --allevents
	
	#if run is successful:
	if  [[ -f "$sample_out/ASGAL.csv" ]];
	then
		#run unification
		echo "Running $tool unificiation..."
		echo "Looking for $tool files in $sample_out/ASGAL/"

		unified_outdir_name="${sample_out}_unmapped_reads_${tool}_dicast_unified"
		echo "Saving unified output to $unified_outdir_name"
		mkdir -p $unified_outdir_name

		anno_file="$workdir/src/ASimulatoR/out/event_annotation.tsv"
		stats_file="${unified_outdir_name}/${fastqname}_unmapped_reads_${tool}_dicast_unified_comparison.txt"

		if [ $combine_events = 0 ];
		then
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -a $sample_out/ASGAL.csv -out $unified_outdir_name -gtf $gtf

			if  [[ -f "$anno_file" ]];
			then
				echo "Running unified comparison..."
				python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${fastqname}_unmapped_reads_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0
			fi
		else
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -a $sample_out/ASGAL.csv -out $unified_outdir_name -gtf $gtf -comb

			if  [[ -f "$anno_file" ]];
			then
				echo "Running unified comparison..."
				python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${fastqname}_unmapped_reads_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0 -comb
			fi
		fi
	else
		echo "Check outputs for ASGAL outputs:"
		echo "$(ls $sample_out)"
		echo "---"
		echo "We should have ASGAL.csv if everything went smooth"
	fi
	
	# if run didn't work out.
	if  [[ ! -f "$sample_out/ASGAL.csv" ]]
	then
		echo "!! ASGAL had a memory issue, recovering from an incomplete run.........." | tee $sample_out/warning.txt
		cat $sample_out/ASGAL/*.events.csv | head -n1 > $sample_out/ASGAL/all.events.csv
		cat $sample_out/ASGAL/*.events.csv | grep -v 'Type,Start' >> $sample_out/ASGAL/all.events.csv
		
		#run unification
		echo "Looking for $tool files in $sample_out/ASGAL/"
		echo "Running $tool unificiation..."

		unified_outdir_name="${sample_out}_unmapped_reads_${tool}_dicast_unified"
		echo "Saving unified output to $unified_outdir_name"
		mkdir -p $unified_outdir_name

		anno_file="$workdir/src/ASimulatoR/out/event_annotation.tsv"
		stats_file="${unified_outdir_name}/${fastqname}_unmapped_reads_${tool}_dicast_unified_comparison.txt"

		if [ $combine_events = 0 ];
		then
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -a $sample_out/ASGAL/all.events.csv -out $unified_outdir_name -gtf $gtf

			if  [[ -f "$anno_file" ]];
			then
				echo "Running unified comparison..."
				python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${fastqname}_unmapped_reads_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0
			fi
		else
			python3 /MOUNT/scripts/unified_output/output_transformer.py create -a $sample_out/ASGAL/all.events.csv -out $unified_outdir_name -gtf $gtf -comb

			if  [[ -f "$anno_file" ]];
			then
				echo "Running unified comparison..."
				python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${fastqname}_unmapped_reads_${tool}_dicast_unified.out -gtf $gtf -stats $stats_file -s -t 0 -comb
			fi
		fi
	fi
done
