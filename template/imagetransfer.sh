#!/bin/bash
echo Saving docker images to the NFS. This may take a while...
for image in $(docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2) 
do echo -ne "\r\033[K Saving the image for $image" && docker save -o /nfs/proj/AS_dockers/images/${image}.tar "proj/${image}:0.01" 
done
echo -ne "\r\033[K All Docker images saved ^_^. Let's begin the export."
echo " "

echo $HOSTNAME
echo logging $USER onto Bert: ; time -p ssh ${USER}@bert 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '
echo logging $USER onto Elmo: ; ssh ${USER}@elmo 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file)  ; done '
echo logging $USER onto Ernie: ; ssh ${USER}@ernie 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '
echo logging $USER onto Grover: ; ssh ${USER}@grover 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '
echo logging $USER onto Kermit: ; ssh ${USER}@kermit 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '

