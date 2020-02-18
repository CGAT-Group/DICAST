#!/bin/bash
second_attempt() {

hisat2 -q  -x /docker_main/latest_human -U "$line"?.fastq -S /myvol1/hisat-output/test.sam --known-splicesite-infile /docker_main/splicesites.txt

}

find /myvol1/ -name "*fastq" -nowarn | sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/hisat-fastqlist

mkdir /myvol1/hisat-output

while read -r line; do



#First attempt
hisat2 -q  -x /docker_main/latest_human -1 "$line"1.fastq -2 "$line"2.fastq -S /myvol1/hisat-output/test.sam --known-splicesite-infile /docker_main/splicesites.txt

trap 'second_attempt $line' ERR

done </myvol1/hisat-fastqlist

wait
rm /myvol1/hisat-fastqlist
echo "script is done"
