#!/bin/bash

for image in $(docker images | grep proj | grep 0.01 | cut -d ' ' -f1 | cut -d '/' -f2) 
do docker save -o /nfs/proj/AS_dockers/images/${image}.tar "proj/${image}:0.01"  
done

echo $HOSTNAME
ssh afenn@bert 'for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done ; echo $HOSTNAME '
ssh afenn@elmo 'for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file)  ; done ; echo $HOSTNAME '
ssh afenn@ernie 'for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done ; echo $HOSTNAME '
ssh afenn@grover 'for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done ; echo $HOSTNAME '
ssh afenn@kermit 'for file in $(find /nfs/proj/AS_dockers/images/ -name '*.tar'); do echo $file $( docker load -i $file) ; done ; echo $HOSTNAME '

