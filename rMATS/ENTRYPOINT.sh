#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=rMATS

mk_outdir $tool
out=${wd}/${output}/${tool:-unspecific}-output  #local variable for readability
tmp_dir=$out/tmp
STAR_index_folder=${wd}/index/star-index/

check_folder_exists(){
	folder=$1
	if [ ! -d "$folder" ];
	then
		echo "Folder $folder didn't exist, creating it..."
		mkdir $folder
		echo "Created $folder"
	else
		echo "Folder $folder already exists."
	fi
}

check_folder_exists $out
check_folder_exists $tmp_dir

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
		echo "Did not find fastq input files in $tmp_dir ; creating them now..."
		create_input_fastq $wd/$casefastq $tmp_dir/s1.txt
		create_input_fastq $wd/$controlfastq $tmp_dir/s2.txt
		echo "Input files created and saved to $tmp_dir."
	elif [[ $use_bam_input_files -eq 1 ]] && ! [[ -f $tmp_dir/b1.txt && -f $tmp_dir/b2.txt ]]
	then
		echo "Did not find bam input files in $tmp_dir ; creating them now..."
		create_input_bam $wd/$casefolder $tmp_dir/b1.txt
		create_input_bam $wd/$controlfolder $tmp_dir/b2.txt
		echo "Input files created and saved to $tmp_dir."
	elif [[ $use_bam_input_files -eq 0 ]] || [[ $use_bam_input_files -eq 1 ]]
	then
		echo "Input files already exist in $tmp_dir."
	else
		echo "Variable use_bam_input_files in the config file is set incorrectly, couldn't recognize value $use_bam_input_files."
		exit 1
	fi

}

rmats_script=$(find /opt/conda/pkgs/rmats*/bin/rmats.py)

if [ ! -d "$STAR_index_folder" ]
then 
	build_STAR_index
	STAR_index_folder=$tmp_dir/STAR_index
else 
	echo STAR-index found, moving on...; fi

create_input_files

if [[ $use_bam_input_files -eq 0 ]]
then
	echo "Start rMATS with fastq input files."
	python $rmats_script --s1 $tmp_dir/s1.txt --s2 $tmp_dir/s2.txt --gtf $wd/$gtf --bi $STAR_index_folder -t paired --readLength $read_length --nthread $ncores --od $out --tmp $tmp_dir
elif [[ $use_bam_input_files -eq 1 ]]
then
	echo "Start rMATS with bam input files."
	python $rmats_script --b1 $tmp_dir/b1.txt --b2 $tmp_dir/b2.txt --gtf $wd/$gtf -t paired --readLength $read_length --nthread $ncores --od $out --tmp $tmp_dir
else
	echo "Variable use_bam_input_files in the config file is set incorrectly, couldn't recognize value $use_bam_input_files."
	exit 1
fi

#cleaner?
