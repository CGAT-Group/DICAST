#!/bin/bash
#Initiating list of images
:> /nfs/proj/AS_dockers/images.txt

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
runscriptenv=0
 ;;
    [nN][oO]|[nN])
 echo "No"
runscriptenv=1
       ;;  
    *)  
 echo "Invalid input..."
 exit 1
 ;;
esac

if [ $runscriptenv -eq 1 ] 
then
newdir=$(zenity --file-selection --directory --text="Where is the ./template/ of the CoMPASS project?")
pushd $newdir
else
pushd ./
fi

# Updating Docker Images
docker-compose build 
echo $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2|wc -l) images to run 
for i in $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2 ); do echo proj/${i}:0.01;done |tee /nfs/proj/AS_dockers/images.txt
echo ------------
echo "Enable dry run? [Y/N] "

read -r -p " Enable dry run? [Y/N] " dryrunner

case $dryrunner in
    [yY][eE][sS]|[yY])
 echo "Yes"
rundry=0
 ;;
    [nN][oO]|[nN])
 echo "No"
rundry=1
       ;;  
    *)  
 echo "Invalid input..."
 exit 1
 ;;
esac

if [ $rundry -eq 0 ] 
then
popd
exit
fi

#Is the runscript executed on the server with SLURM? 
if command -v squeue &> /dev/null #checking if SLURM exists.
then
read -r -p "Images will not be transfered, if the clusters are updated. Are the clusters updated? [Y/N] " input

case $input in
    [yY][eE][sS]|[yY])
 echo "Yes"
updatecluster=0
 ;;
    [nN][oO]|[nN])
 echo "No"
updatecluster=1
       ;;
    *)
 echo "Invalid input..."
 exit 1
 ;;
esac
if [ $updatecluster -eq 1 ]
then
/bin/bash ./imagetransfer.sh
fi
sbatch ./slurm-docker-runner.sh
else
#Runscript without SLURM
bash ./docker-runner.sh
fi
popd
