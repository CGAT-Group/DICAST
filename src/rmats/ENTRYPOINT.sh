#!/bin/bash

tool=rmats
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

mk_outdir
tmp_dir=tmp

create_input_fastq(){
	folder=$1
	output_file=$2
	find $folder/ -maxdepth 1 -type f -name '*.fastq' | awk -F $fastqpair1suffix"|"$fastqpair2suffix '{print $1}' | sort |
		uniq | awk '{new_str=$1suffix1":" $1suffix2","; printf new_str}' suffix1=$fastqpair1suffix suffix2=$fastqpair2suffix  >> $output_file
	sed -i '$s/,$//' $output_file
}

create_input_bam(){
	folder=$1
	output_file=$2
	find $folder/ -maxdepth 1 -type f -name '*.bam' | awk -v ORS=, '{ print $0 }' | sed '$s/,$//' >> $output_file
}

create_input_files(){
	if [[ $use_bam_input_files -eq 0 ]] && ! [[ -f $tmp_dir/s1.txt && -f $tmp_dir/s2.txt ]]
	then
		create_input_fastq $casefastq $tmp_dir/s1.txt
		create_input_fastq $controlfastq $tmp_dir/s2.txt
	elif [[ $use_bam_input_files -eq 1 ]] && ! [[ -f $tmp_dir/b1.txt && -f $tmp_dir/b2.txt ]]
	then
		create_input_bam $casefolder $tmp_dir/b1.txt
		create_input_bam $controlfolder $tmp_dir/b2.txt
	elif [[ $use_bam_input_files -eq 0 ]] || [[ $use_bam_input_files -eq 1 ]]
	then
		echo "Input files already exist in $tmp_dir."
	else
		echo "Variable use_bam_input_files in the config file is set incorrectly, couldn't recognize value $use_bam_input_files."
		exit 1
	fi

}

rmats_script=$(find /opt/conda/pkgs/rmats*/bin/rmats.py)
check_star_index
create_input_files

if [[ $use_bam_input_files -eq 0 ]]
then
	echo "Start rMATS with fastq input files."
	python $rmats_script --s1 $tmp_dir/s1.txt --s2 $tmp_dir/s2.txt --gtf $gtf --bi $star_index -t paired --readLength $read_length --nthread $ncores --od $outdir --tmp $tmp_dir
elif [[ $use_bam_input_files -eq 1 ]]
then
	echo "Start rMATS with bam input files."
	python $rmats_script --b1 $tmp_dir/b1.txt --b2 $tmp_dir/b2.txt --gtf $gtf -t paired --readLength $read_length --nthread $ncores --od $outdir --tmp $tmp_dir
else
	echo "Variable use_bam_input_files in the config file is set incorrectly, couldn't recognize value $use_bam_input_files."
	exit 1
fi

#cleaner?
