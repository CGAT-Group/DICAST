#!/bin/bash

tool=asimulator

config_path=$(pwd)/scripts          #when executed by Snakefile, $pwd = working directory of snakemake.
echo "Config folder path: $config_path."
asimulator_inputdir=$(pwd)/src/ASimulatoR/in
echo "ASimulatoR input path: $asimulator_inputdir."
asimulator_outputdir=$(pwd)/src/ASimulatoR/out
echo "ASimulatoR output path: $asimulator_outputdir."

# acknowledging config set for DICAST
source $config_path/config.sh
source $config_path/mapping_config.sh

# prepping input dir for ASimulator @ src/ASimulator/in.
bowtie_fastadirname=$(echo $bowtie_fastadir | sed 's/\/MOUNT\///g')
for i in $(ls $bowtie_fastadirname); do ln $bowtie_fastadirname/$i $asimulator_inputdir/$i ; done
ln $( echo $fasta | sed 's/\/MOUNT\///g') $asimulator_inputdir/$(basename $fasta)
ln $( echo $inputdir| sed 's/\/MOUNT\///g')/$asimulator_gtf $asimulator_inputdir/$asimulator_gtf
cp $config_path/ASimulatoR_config.R $asimulator_inputdir/runASS.R

#clearing output dir for ASimulator @ src/ASimulator/out
echo Removing old ASimulatoR runs from $asimulator_outputdir
rm $asimulator_outputdir/* -R

# checking if ASimulatoR docker exists.
if [[ "$(docker inspect --type=image biomedbigdata/asimulator 2>/dev/null)" == "[]" ]]; then
  echo "Pulling ASimulatoR docker image ..."
  docker pull biomedbigdata/asimulator
else
  echo "Found ASimulatoR docker image locally."
fi

docker run --rm --name dicast-$tool --user $(id -u):$(id -g) -v $asimulator_inputdir:/input -v $asimulator_outputdir:/output biomedbigdata/asimulator

# Bringing the outputs to inputdir.

ln $asimulator_outputdir/*/*gtf $( echo $inputdir| sed 's/\/MOUNT\///g')/ASimulatoR.gtf
ln $asimulator_outputdir/*/*fastq $( echo $fastqdir| sed 's/\/MOUNT\///g')
ln $asimulator_outputdir/*/*gff3 $( echo $inputdir| sed 's/\/MOUNT\///g')/ASimulatoR.gff3