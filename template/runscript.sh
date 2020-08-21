#!/bin/bash
:> /nfs/proj/AS_dockers/images.txt
#Warning
{
i=1
sp="/-\|"

printf "\
=========================================================================================\n\
=  Welcome to CoMPASS  built by the Alternative Splicing - Docker Team at TUM =\n\
=========================================================================================\n\
\n\
If you need help with this Pipeline, please contact: \n\
Amit FENN, Chair of Experimental Bioinformatics, Technical University of Munich (TUM), Freising, Germany \n\
\n\
The Docker team: Fanny Rößler, Alexander Dietrichi, Tim Faro and Amit\
\n\
If you liked our Pipeline, then consider buying us a beer.\n \n \n"

echo "This script requires to be executed in the folder ./template/ of the CoMPASS project"
while [ $i -le 70 ]
do
printf "\b${sp:i++%${#sp}:1}"
sleep 0.2
done
}
docker-compose build 
echo $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2|wc -l) images to run 
for i in $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2 ); do echo proj/${i}:0.01;done |tee /nfs/proj/AS_dockers/images.txt
echo ------------
echo " "
read -r -p "update clusters? [Y/n] " input

case $input in
    [yY][eE][sS]|[yY])
 echo "Yes"
updatecluster=1
 ;;
    [nN][oO]|[nN])
 echo "No"
updatecluster=0
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
