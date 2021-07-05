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
if [ ! -z "$asimulator_outputdir" ]
    then rm $asimulator_outputdir/* -R
    else echo asimulator_outputdir variable unset. Please correct src/run_asimulator.sh . && exit 1
fi

# checking if ASimulatoR docker exists.
if [[ "$(docker inspect --type=image biomedbigdata/asimulator 2>/dev/null)" == "[]" ]]; then
  echo "Pulling ASimulatoR docker image ..."
  docker pull biomedbigdata/asimulator
else
  echo "Found ASimulatoR docker image locally."
fi

docker run --rm --name $USER-$RANDOM-dicast-$tool --user $(id -u):$(id -g) -v $asimulator_inputdir:/input -v $asimulator_outputdir:/output biomedbigdata/asimulator

# Bringing the outputs to inputdir.
set +e
mv $asimulator_outputdir/*gtf $( echo $inputdir| sed 's/\/MOUNT\///g')/ASimulatoR.gtf
mv $asimulator_outputdir/*gff3 $( echo $inputdir| sed 's/\/MOUNT\///g')/ASimulatoR.gff3
for i in $(ls ${asimulator_outputdir}/*fastq)
  do mv ${asimulator_outputdir}/$i $( echo $fastqdir| sed 's/\/MOUNT\///g')/Simulated-$( basename $i )
  done