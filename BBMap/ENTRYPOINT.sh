#!/bin/bash
second_attempt() {

#Second attempt isn't working for BBMap. :(
/bbmap/bbmap.sh ref=$(find /docker_main/ -name "*Homo*.fa*") 

/bbmap/bbmap.sh in="$line"?.fastq outm=/myvol1/BBmap-output/"${line##*/}"_mapped.sam outu=/myvol1/BBmap-output/"${line##*/}"_unmapped.sam
}

#START here: Make a list of fastq files
find /myvol1/ -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/BBmap-fastqlist

#make output directories
mkdir -p /myvol1/BBmap-output/
cd /myvol1/BBmap-output/

#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
/bbmap/bbmap.sh ref=$(find /docker_main/ -name "*Homo*.fa*")

/bbmap/bbmap.sh in1="$line"1.fastq in2="$line"2.fastq outm=/myvol1/BBmap-output/"${line##*/}"_mapped.sam outu=/myvol1/BBmap-output/"${line##*/}"_unmapped.sam

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/BBmap-fastqlist

#wait for all processes to end
wait

#cleaning up
rm /myvol1/BBmap-fastqlist
echo "script is done"

