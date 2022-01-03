#!/bin/bash
tools=( asgal aspli eventpointer irfinder majiq sgseq spladder whippet subjunc bbmap contextmap crac dart gsnap hisat minimap mapsplice segemehl star )
spin='-\|/'

dockerhub_prefix="dicastproj/dicast"
locally_built_prefix="dicast/"
echo Uninstalling DICAST

for image in ${tools[*]}; do  docker rmi $dockerhub_prefix:$image 2> /dev/null ; docker rmi ${locally_built_prefix}${image}:0.2 2> /dev/null ; done &
pid=$! # Process Id of the previous running command
i=0
while kill -0 $pid 2>/dev/null
do
  i=$(( (i+1) %4 ))
  printf "\r${spin:$i:1}"
  sleep .1
done

conda activate base 2> /dev/null
conda remove --name dicast-snakemake --all -y &
pid=$! # Process Id of the previous running command
i=0
while kill -0 $pid 2>/dev/null
do
  i=$(( (i+1) %4 ))
  printf "\r${spin:$i:1}"
  sleep .1
done
echo " You can now delete this git folder and everything it contains. Consider running a docker images prune"