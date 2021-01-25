#!/bin/bash


#read update_tools.txt
file="update_tools.txt"
branch="develop"

while read tool
do	
	echo "Building image for $tool"
	cd ../src/$tool
	#change tool-name to lowercase
	tool_lc=$(echo $tool | tr '[:upper:]' '[:lower:]')
	echo $tool_lc
	#build docker-image for this tool
	docker build -t gitlab.lrz.de:5005/ge46ban/dockers/$branch/$tool_lc:latest .
	#push image into registry
	docker push gitlab.lrz.de:5005/ge46ban/dockers/$branch/$tool_lc:latest
	cd ../../scripts
	echo "----------------------------"
	echo "-----------NEXT-------------"
	echo "----------------------------"
done <$file
