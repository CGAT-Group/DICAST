#!/bin/bash
#Initiating list of images
if command -v squeue &> /dev/null #checking if SLURM exists.
then
:> /nfs/proj/AS_dockers/images.txt
SLURM=1
else
:> ./dockerrunlist
SLURM=0
fi

#Welcome message
printf "\
=========================================================================================\n\
=  Welcome to CoMPASS  built by the Alternative Splicing - Docker Team at TUM =\n\
=========================================================================================\n\
\n\
If you need help with this Pipeline, please contact: \n\
Amit FENN, Chair of Experimental Bioinformatics, Technical University of Munich (TUM), Freising, Germany \n\
\n\
The Docker team: Fanny Rößler, Alexander Dietrich, Tim Faro and Amit\
\n\
If you liked our Pipeline, then consider buying us a beer.\n \n \n"
echo "This script requires to be executed in the folder ./template/ of the CoMPASS project"

# Checking Run Context
read -r -p "Is this the ./template/ folder of the CoMPASS project [Y/N]? " confirmdir
case $confirmdir in
[yY][eE][sS]|[yY])
echo "Yes"
pushd ./
;;
[nN][oO]|[nN])
echo "No"
newdir=$(zenity --file-selection --directory --text="Where is the ./template/ of the CoMPASS project?")
pushd $newdir
;;  
*)  
echo "Invalid input..."
exit 1
;;
esac

# Updating Docker Images
docker-compose build 
#echo $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2|wc -l) images to run 
#Is the runscript executed on the server with SLURM? 
if command -v squeue &> /dev/null #checking if SLURM exists.
then
read -r -p "Images will not be transfered, if the clusters are updated. Are the clusters updated? [Y/N] " input

case $input in
[yY][eE][sS]|[yY])
echo "Yes"
;;
*)
echo "No"
/bin/bash ./imagetransfer.sh
;;
esac

read -r -p " Enable dry run? [Y/N] " dryrunner
case $dryrunner in
[yY][eE][sS]|[yY])
echo "Yes"
exit 1
popd
;;
*)  
echo "No"
;;  
esac

# Select Dockers to run
IFS='|' read -ra dockerarray <<<  $(zenity --list --height 800 --width 400  --checklist --title="Select dockers to run"\
    --text="Select the tools to run"\
    --column="Use"\
    --column="Docker"\
    $(for i in $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2 ); do echo  TRUE proj/${i}:0.01;done ))
unset $IFS
case $SLURM in
1)	#Print the split string
for i in "${dockerarray[@]}"
do
    echo ${i} | tee -a  /nfs/proj/AS_dockers/images.txt 
done
;;
0)  #Print the split string 
for i in "${dockerarray[@]}"
do
    echo ${i} | tee -a  ./dockerrunlist
done
;;
esac
#for i in $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2 ); do echo proj/${i}:0.01;done |tee /nfs/proj/AS_dockers/images.txt ./dockerrunlist
# Allow Exit
echo ------------


if command -v squeue &> /dev/null #checking if SLURM exists.
then
read -r -p "Run Sbatch? " sinput

case $sinput in
[yY][eE][sS]|[yY])
echo "Yes"
bash arraymaker.sh
sbatch ./slurmsubmit.sh
;;
*)
echo "No"
exit 1
;;
esac
fi
else
echo Running dockers locally.
#Runscript without SLURM
bash ./docker-runner.sh
fi
popd
