#!/bin/bash

#Read SAM and BAM files, first within the path described in a parameter, else in $samfolder/$bamfolder

readsamfiles(){
	find $wd${1:-$samfolder} -name "*.sam" -nowarn > $wd/tmp/$tool-samlist 
		chmod  777 $wd/tmp/$tool-samlist 
}

readbamfiles(){
	find $wd${1:-$bamfolder} -name "*.bam" -nowarn > $wd/tmp/$tool-bamlist
		chmod  777 $wd/tmp/$tool-bamlist
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




#cleaning up 
cleaner(){
	rm $wd/tmp/$tool-?amlist;
	echo "script is done"; 
	chmod -R 777 $wd/$out 
} 
