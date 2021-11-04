#!/bin/bash

#Read gtf file and test if its available
#Parameter: Full path to gtf file

test_gtf(){
	if [[ ! -f $1 ]]
	then
		echo File not found. Check the path for the ${gtf} gtf file.
		exit 1
	else
		echo found gtf file, moving on...
	fi
}


#Read fasta file and test if its available
#Parameter: Full path to fasta file

test_fasta(){
        if [[ ! -f $1 ]]
        then
                echo File not found. Check the path for the ${fasta} fasta file.
                exit 1
        else
                echo found fasta file, moving on...
        fi
}




#test if BAM file is available and has index file with it
#Parameter: Full path to BAM-file
test_bam(){
	if [[ ! -f $1  ]]
	then
		echo "File not found. Check path for the bam file(s): is it in ${bamdir}?"
		exit 1
	else
		echo found bam file, moving on...
		if [[  ! -f ${1}.bai  ]]
		then
			echo File not found. Did not find index file for ${1}. Check if it has the same name as ${1} and is in ${bamdir}.
		else
			echo found bam index, moving on...
		fi
	fi
}


#take all BAM & BAM-index files from bamfolder(or case-/controlfolder) and store paths in the new bamlist contained in the tool-output folder
#Parameter1: folder to check for bam files;
#Parameter2: optional output filename (standard is bamlist)
readbamfiles(){
	#if bamlist exists: delete
	rm /tmp/${2:-bamlist} 2>/dev/null
	find ${1:-$controlbam} -maxdepth 1 -name "*.bam" -nowarn -not -empty >> /tmp/${2:-bamlist}
	chmod  777 /tmp/${2:-bamlist}
}

#just as readbamfiles, but for SAM-files
#parameter1: optional different bamdirectory
readsamfiles(){
	#touch $wd/$output/${tool:-unspecific}-output/samlist
	find ${1:-$controlbam} -maxdepth 1 -name "*.sam" -nowarn -not -empty >> /tmp/samlist
	chmod  777 /tmp/samlist
}



#sort bams if unsorted
sortnindexbam(){
	if [ ! -s ${1}.bai ]
	#if [[ ! "$(samtools view -H $1 | grep SO: | cut -f 3 | cut -d ":" -f2)" == "coordinate" ]]; then
		#local bamname=$(basename -s .bam $1)
		then	samtools sort $1 -o $1 && \
			samtools index $1  -@ 4 ;
	fi
}


#Function: Check if bam exists in sam folder or bam folder, if not, build a bam from sam in parameter, sort and index it.
#Parameter1: Full file path to sam file
#Parameter2: directory on which this function is used ($controlbam by default or $casebam & $controlbam); default is $controlbam
makebamfromsam(){
	local sampath=$(dirname $1)
	local samfileprefix=$(basename -s .sam $1)

		if [[ ! -e ${sampath}/${samfileprefix}.bam ]]
			then
				if [[ ! -e ${2:-$controlbam}/${samfileprefix}.bam ]]
					then
					{
#Make a Bam file
						echo making bam of $1, in ${2:-$controlbam}/$samfileprefix.bam.  This may take a while..
							samtools view -b -S -T $fasta $1 > ${2:-$controlbam}/${samfileprefix}.bam
#Sort bam file
						echo sorting and indexing bam file $samfileprefix.bam
							sortnindexbam "${2:-$controlbam}/$samfileprefix.bam"
					}
				else
					echo bam file for $1 exists in ${2:-$controlbam}
						sortnindexbam "${2:-$controlbam}/$samfileprefix.bam"
				fi
		else
			echo bam file for $1 exists in $sampath
				#if [[ ! -e ${2:-$controlbam}/${samfileprefix}.bam ]] ; then
				#	mv $1 "${2:-$controlbam}"
				#		echo "moved $1 to :" ${2:-$controlbam}
				#		sortnindexbam "${2:-$controlbam}/${samfileprefix}.bam"
				#fi
		fi
#readbamfiles
readbamfiles $controlbam controlbamlist
readbamfiles $casebam casebamlist
}

#function to handle sam files in either bamdir (for as_tools) or case/control-bamdir (for ds_tools)
#parameter: 0 to use for as_tools; 1 to use in ds_tools
handlesamfiles(){
	#clear samlist file first
	rm -f /tmp/samlist 2>/dev/null
	rm -f /tmp/bamlist 2>/dev/null
	if [[ $1 = 0 ]]
	then
		echo "Looking for SAM files in $controlbam and converting them to BAM-files..."
		readsamfiles $controlbam
		readbamfiles $controlbam
		#make bam file out of all samfiles in samlist
		for filename in $(cat /tmp/samlist)
		do
		        makebamfromsam $filename $controlbam
		done

	else
		echo "Looking for SAM files in $casebam and $controlbam and converting them to BAM-files..."
		readsamfiles $casebam
		readbamfiles $casebam
		#make bam file out of all samfiles in samlist
                for filename in $(cat /tmp/samlist)
                do
                        makebamfromsam $filename $casebam
                done
		echo "-------------------"
		#clear samlist again
		rm -f /tmp/samlist 2>/dev/null
		readsamfiles $controlbam
                #make bam file out of all samfiles in samlist
                for filename in $(cat /tmp/samlist)
                do
                        makebamfromsam $filename $controlbam
                done

	fi

}



# read fastq files in fastqfolder and save paths ins fastqlist file
readfastqs(){
        find ${1:-$fastqdir} -maxdepth 2 -name "*.fastq" -nowarn > /tmp/fastqlist
        chmod  777 /tmp/fastqlist
		while read -r line; do
			#mkdir -p $outdir/$(basename $(dirname $line))/
			chmod -R 777 $outdir/
		done </tmp/$tool-fastqlist
}


#make the output directory: $output/$tool-output
#Parameter: the tools name
 mk_outdir(){
 	mkdir -p $outdir/
 	chmod -R 777 $outdir
 }


#make the output folder for the specific sample (BAM-file) and return the folder path
#Parameter: BAM filename of sample
mk_sample_out(){
	tmp="${1##*/}"	#get basename of file
	sample_out="${tmp%%.*}"	#remove all file extensions after first.
	if [ "$tool" = "asgal" ] || [ "$tool" = "irfinder" ]
		then 
		mkdir -p $outdir/${sample_out}_unmapped
		chmod -R 777 $outdir/${sample_out}_unmapped
		echo $outdir/${sample_out}_unmapped
		 else
		mkdir -p $outdir/${sample_out}_output
		chmod -R 777 $outdir/${sample_out}_output
		echo $outdir/${sample_out}_output
	fi
}


build_STAR_index(){
	mkdir -p $star_index
	STAR --runMode genomeGenerate --genomeDir $star_index --genomeFastaFiles $fasta --sjdbGTFfile $gtf --runThreadN $ncores --outFileNamePrefix $star_index/Star_mapped_ --sjdbOverhang 100
	wait
	echo "Built STAR index and saved into $star_index ..."
}

check_star_index(){
	echo "Looking for STAR index in $star_index ..."
	if [ ! -d $star_index ]; then
		echo "Did not find index, building it now into $star_index ..."
		build_STAR_index
	else
		echo "Found STAR index, moving on ..."
	fi
}


#cleaning up
cleaner(){
	rm -f $outdir/bamlist 2>/dev/null
	rm -f $outdir/samlist 2>/dev/null
	rm -f $outdir/fastqlist 2>/dev/null
	rm -rf $outdir/tmp 2>/dev/null
	echo script is done
	chmod 777 $outdir -R
	#exit
}






###### functions used for tools which can also do differential splicing analysis

#combine files of case & control folders into one tmp/case_control folder and return the path
#Parameter: 1: path of casefolder 2: path of controlfolder
combine_case_control(){
	mkdir -p $outdir/case_control
	chmod -R 777 $outdir/case_control

	cp -R $1/. $outdir/case_control
	cp -R $2/. $outdir/case_control

	echo $outdir/case_control
}


cleaner_diff(){
	rm -rf $outdir/case_control 2>/dev/null
	rm -f $outdir/bamlist 2>/dev/null
       	echo script is done.
	exit
}

#################
#    Logging    #
#################

echo_vars() {
	echo "############################################"
	echo "VARIABLE SETTINGS"
	echo -e """

tool:\t\t$tool

DIRECTORIES
workdir:\t$workdir
outdir: \t$outdir
inputdir:\t$inputdir
indexdir:\t$indexdir
fastqdir:\t$fastqdir	# only for asgal, irfinder, kissplice, whippet
star_index:\t$star_index	# only for kissplice

FILES
fasta:\t\t$fasta	# only for asgal, irfinder, whippet
gtf:\t\t$gtf	# not for majiq
gff:\t$gff	# only for majiq

PARAMETERS
ncores: \t$ncores	# only for applicable tools
read_length:\t$read_length	# only for applicable as_tools
use_bam_input_files:\t$use_bam_input_files	# only for irfinder
differential:\t$differential	# 0 = non differential; 1 = differential

CONTROL (Default for NON DIFFERENTIAL)
controlfolder:\t$controlfolder
controlbam:\t$controlbam
controlfastq:\t$controlfastq
controlprefix:\t$controlprefix	# only for kissplice

CASE
casefolder:\t$casefolder
casebam:\t$casebam
casefastq:\t$casefastq
caseprefix:\t$caseprefix	# only for kissplice
"""
	echo "############################################"
	echo
}

start_logging() {
	mkdir -p $outdir/logs
	current_time=$(date "+%Y.%m.%d_%H:%M:%S")
	log_file=$outdir/logs/${tool}_${current_time}.log
	touch $log_file
	echo -e "\nlogs will be stored in ${log_file}\n"
	exec &> >(tee -a -i "$log_file")
	echo_vars
	#exec 2> >(tee -a -i "${log_file}")
	#exec >> "${log_file}"

}
