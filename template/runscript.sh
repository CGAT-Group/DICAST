#!/bin/bash
:> /nfs/proj/AS_dockers/images.txt
echo "This script requires to be executed in the folder ./template/ of the CoMPASS project" ;echo -ne '.\r';sleep 1;echo -ne '..\r';sleep 1; echo -ne '...\r'; echo -ne '\n';sleep 1
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
