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
							samtools view -bS $1 > ${bamfolder}/${samfileprefix}.bam
					} && {
#Sort bam file
						echo sorting bam file $samfileprefix.bam
							samtools sort "${bamfolder}/$samfileprefix.bam" -o "${bamfolder}/$samfileprefix-sorted.bam" && \
							rm "$wd/${bamfolder}/$samfileprefix.bam"

					} &&{
#Index bam file
						echo indexing bam file $samfileprefex-sorted.bam
							samtools index "$samfileprefix-sorted.bam" -@ 4
					}
				else 
					echo bam file for $1 exists in $bamfolder
						fi
		else
			echo bam file for $1 exists in $sampath
				fi
} 




#cleaning up 
cleaner(){
	rm $wd/tmp/$tool-?amlist;
	echo "script is done"; 
	chmod -R 777 $wd/$out 
} 
