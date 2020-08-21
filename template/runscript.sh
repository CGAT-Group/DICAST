#!/bin/bash
:> /nfs/proj/AS_dockers/images.txt
for i in $( docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2 ); do echo proj/${i}:0.01;done |tee /nfs/proj/AS_dockers/images.txt

read -r -p "update clusters? [Y/n] " input

case $input in
    [yY][eE][sS]|[yY])
 echo "Yes"
updatecluster=1
 ;;
    [nN][oO]|[nN])
 echo "No"
       ;;
    *)
 echo "Invalid input..."
 exit 1
 ;;
esac

if [ $updatecluster -eq 1 ]
then
/bin/bash /nfs/home/users/afenn/Dockers/dockers/template/imagetransfer.sh
fi
sbatch /nfs/home/users/afenn/Dockers/dockers/template/slurm-docker-runner.sh
