#!/bin/bash 

#Read gtf file and test if its available
#Parameter: Full path to gtf file

test_gtf(){
	if [[ ! -f $1 ]]
	then 
		echo File not found. Check the path for the ${gtf} gtf file: is it under ${wd}/${gtf}?
		exit 1
	else
		echo found gtf file, moving on...
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
 

#TODO: functions for fastq-files (list them + check for correct path)
#test_fastq(){}
#list_fastqs(){}


#make the output directory: $output/$tool-output
#Parameter: the tools name
mk_outdir(){
	mkdir -p /$wd/$output/${1:-unspecific}-output
	chmod -R 777 /$wd/$output/${1:-unspecific}-output
}

#cleaning up
cleaner(){
	echo script is done
	exit
}






###### functions used for tools which can also do differential splicing analysis

#combine files of case & control folders into one tmp/case_control folder and return the path
combine_case_control(){
	mkdir -p $wd/tmp/case_control
	chmod -R 777 $wd/tmp/case_control

	cp -r $wd/casefolder/ $wd/tmp/case_control/
	cp -r $wd/controlfolder/ $wd/tmp/case_control/

	echo ${wd}/tmp/case_control
}


cleaner_diff(){
	rm -rf $wd/tmp/case_control
       	echo script is done.
	exit
}
