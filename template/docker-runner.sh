#!/bin/bash
while read image
do;
{docker run -v "/nfs/proj/AS_dockers/:/myvol1" --user $(id -u):$(id -g) --rm $image ; echo This docker exited with status is: $?}> ./runlog-${image}.log
echo $image $USER
done < ./dockerrunlist 
