#!/bin/bash

tools=( sgseq irfinder )
#tools=( eventpointer asgal )
#tools=(  mapsplice segemehl star aspli irfinder sgseq spladder whippet asgal eventpointer majiq subjunc minimap hisat gsnap dart crac contextmap bbmap)
#tools=( star  aspli  irfinder  sgseq spladder whippet asgal eventpointer majiq )
#unused tools: subjunc bbmap contextmap crac dart gsnap hisat minimap mapsplice segemehl asgal eventpointer majiq
echo tool-list currently: ${tools[*]}

gitlabimagename="gitlab.lrz.de:5005/ge46ban/dockers/develop/"
dockerhubimagename="dicastproj/dicast"


dockerlist=$(docker images | grep dicast/ | cut -d ' ' -f1)

for i in ${tools[*]}
do echo tool: $i
docker tag dicast/$i:0.2 ${gitlabimagename}${i}:0.3
docker push ${gitlabimagename}${i}:0.3
done

#docker-compose -f /nfs/proj/AS_dockers/amit/dockers/scripts/snakemake/docker-compose.yml build

for i in ${tools[*]}
do echo tool: $i
docker tag dicast/$i:0.2 ${dockerhubimagename}:${i}
docker push ${dockerhubimagename}:${i}
done


#docker-compose -f /nfs/proj/AS_dockers/amit/dockers/scripts/snakemake/docker-compose.yml build

echo Moving on to delete images
docker images | grep $gitlabimagename
docker images | grep $dockerhubimagename

#	echo Countdown ends at 20.
#for i in $(seq 1 2)
#do sleep 1
#echo $i
#done

	for i in ${tools[*]}
	do echo tool: $i
	echo docker rmi ${gitlabimagename}${i}:0.3
	done

	for i in ${tools[*]}
	do echo tool: $i
	echo docker rmi ${dockerhubimagename}:${i}
	done
