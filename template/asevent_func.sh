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
		echo File not found. Check path for the $bamfile file: is it in ${wd}/${bamfolder}?
		exit 1
	else
		echo found bam file, moving on...
		if [[  ! -f ${1}.bai  ]]
		then
			echo File not found. Did not find index file for ${1}. Check if it has the same name as ${1} and is in ${wd}/${bamfolder}.
		else
			echo found bam index, moving on...
		fi
	fi
}


#take all BAM & BAM-index files from bamfolder(or case-/controlfolder) and store paths in the new bamlist contained in the tool-output folder
#Parameter: folder to check for bam files; 
readbamfiles(){
	find $wd/${1:-$bamfolder} -maxdepth 1 -name "*.bam" -nowarn > $wd/$output/${tool:-unspecific}-output/bamlist
	chmod  777 $wd/$output/${tool:-unspecific}-output/bamlist
}

#just as readbamfiles, but for SAM-files
readsamfiles(){
	#touch $wd/$output/${tool:-unspecific}-output/samlist
	find $wd/${1:-$bamfolder} -maxdepth 1 -name "*.sam" -nowarn > $wd/$output/${tool:-unspecific}-output/samlist
	chmod  777 $wd/$output/${tool:-unspecific}-output/samlist
}


indexbam(){
	if [ ! -s ${1}.bai ] ; then 
		samtools index $1  -@ 4
	fi
}

#sort bams if unsorted
sortnindexbam(){
	if [[ ! "$(samtools view -H $1 | grep SO: | cut -f 3 | cut -d ":" -f2)" == "coordinate" ]]; then
		local bamname=$(basename -s .bam $1)
			samtools sort $1 -o "${bamname}.bam" && \
			indexbam "${bamname}.bam"
	fi
}


#Function: Check if bam exists in sam folder or bam folder, if not, build a bam from sam in parameter, sort and index it.
#Parameter: Full file path to sam file
#
makebamfromsam(){
	local sampath=$(dirname $1)
	local samfileprefix=$(basename -s .sam $1)

		if [[ ! -e ${sampath}/${samfileprefix}.bam ]] 
			then
				if [[ ! -e ${wd}/${bamfolder}/${samfileprefix}.bam ]]
					then
					{
#Make a Bam file
						echo making bam of $1, in $wd/$bamfolder/$samfileprefix.bam.  This may take a while..
							samtools view -bS $1 > $wd/${bamfolder}/${samfileprefix}.bam
#Sort bam file
						echo sorting and indexing bam file $samfileprefix.bam
							sortnindexbam "$wd/${bamfolder}/$samfileprefix.bam"
					}
				else 
					echo bam file for $1 exists in $bamfolder
						sortnindexbam "$wd/${bamfolder}/$samfileprefix.bam"				
				fi
		else
			echo bam file for $1 exists in $sampath
				if [[ ! -e ${wd}/${bamfolder}/${samfileprefix}.bam ]] ; then
					mv $1 "$wd/${bamfolder}" 
						echo "moved $1 to :" $bamfolder
						sortnindexbam "${wd}/${bamfolder}/${samfileprefix}.bam"
				fi
		fi
} 


# read fastq files in fastqfolder and save paths ins fastqlist file
readfastqs(){
        find $wd/${1:-$fastqfolder} -maxdepth 1 -name "*.fastq" -nowarn > $wd/$output/${tool:-unspecific}-output/fastqlist
        chmod  777 $wd/$output/${tool:-unspecific}-output/fastqlist
}


#make the output directory: $output/$tool-output
#Parameter: the tools name
mk_outdir(){
	mkdir -p $wd/$output/${1:-unspecific}-output
	chmod -R 777 $wd/$output/${1:-unspecific}-output
}


#make the output folder for the specific sample (BAM-file) and return the folder path
#Parameter: BAM filename of sample
mk_sample_out(){
	tmp="${1##*/}"	#get basename of file
        sample_out="${tmp%%.*}"	#remove all file extensions after first . 
        mkdir -p $wd/$output/${tool:-unspecific}-output/$sample_out-output
	chmod -R 777 $wd/$output/${tool:-unspecific}-output/$sample_out-output

	echo $wd/$output/${tool:-unspecific}-output/$sample_out-output
}

#cleaning up
cleaner(){
	rm -f $wd/$output/${tool:-unspecific}-output/bamlist
	rm -f $wd/$output/${tool:-unspecific}-output/samlist
	rm -f $wd/$output/${tool:-unspecific}-output/fastqlist
	rm -rf $wd/$output/${tool:-unspecific}-output/tmp
	echo script is done
	#exit
}






###### functions used for tools which can also do differential splicing analysis

#combine files of case & control folders into one tmp/case_control folder and return the path
#Parameter: 1: path of casefolder 2: path of controlfolder
combine_case_control(){
	mkdir -p $wd/$output/${tool:-unspecific}-output/case_control
	chmod -R 777 $wd/$output/${tool:-unspecific}-output/case_control

	cp -R $wd/$1/. $wd/$output/${tool:-unspecific}-output/case_control
	cp -R $wd/$2/. $wd/$output/${tool:-unspecific}-output/case_control

	echo $wd/$output/${tool:-unspecific}-output/case_control
}


cleaner_diff(){
	rm -rf $wd/$output/${tool:-unspecific}-output/case_control
	rm -f $wd/$output/${tool:-unspecific}-output/bamlist
       	echo script is done.
	exit
}
