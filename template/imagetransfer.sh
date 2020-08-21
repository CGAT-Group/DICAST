#!/bin/bash

for image in $(docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2) 
do docker save -o /nfs/proj/AS_dockers/images/${image}.tar "proj/${image}:0.01"  
done

echo $HOSTNAME
echo logging onto Bert: ; ssh afenn@bert 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '
echo logging onto Elmo: ; ssh afenn@elmo 'echo $HOSTNAME: for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file)  ; done '
echo logging onto Ernie: ; ssh afenn@ernie 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '
echo logging onto Grover: ; ssh afenn@grover 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '
echo logging onto Kermit: ; ssh afenn@kermit 'echo $HOSTNAME; for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done '

