#!/bin/bash

tool=asimulator

print_help(){
  echo "$tool run script"
  echo " "
  echo "Usage: run_asimulator.sh [OPTIONS]"
  echo " "
  echo "Options:"
  echo "-h    show this help message and exit."
  echo "-c    specify the path to the config folder containing the config.sh and ASimulatoR_config.R file."
  exit 0
}

if [[ "$#" -eq "0" ]]; then
   print_help
elif [[ "${1:0:1}" != "-" ]]; then
  echo "Argument needs flag! (Missing hyphen?)"
  echo " "
  print_help
fi

config_path=""
while getopts ":hc:" opt; do
  case ${opt} in
    h )
      print_help
      ;;
    c )
      config_path=$OPTARG
      shift
      ;;
    \? )
      echo "Couldn't parse argument $OPTARG!"
      echo " "
      print_help
      ;;
  esac
done

echo "Config folder path: $config_path."

source $config_path/config.sh
echo $asimulator_inputdir
echo $asimulator_outputdir

cp $config_path/ASimulatoR_config.R $asimulator_inputdir/runASS.R

echo "Pulling ASimulatoR docker image ..."
docker pull biomedbigdata/asimulator

docker run --user $(id -u):$(id -g) -v $asimulator_inputdir:/input -v $asimulator_outputdir:/output biomedbigdata/asimulator
