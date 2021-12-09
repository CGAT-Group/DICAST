#!/bin/bash
tools=( asgal aspli eventpointer irfinder majiq sgseq spladder whippet subjunc bbmap contextmap crac dart gsnap hisat minimap mapsplice segemehl )

dockerhub_prefix = "dicastproj/dicast"
locally_built_prefix = "dicast/"

for image in ${tools[*]}; do  docker rmi $dockerhub_prefix:$image && docker rmi ${locally_built_prefix}${image}:0.2 ; done
conda activate base
conda remove --name dicast-snakemake --all
echo " You can now delete this git folder and everything it contains."