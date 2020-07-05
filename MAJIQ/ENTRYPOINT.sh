#!/bin/bash

source /myvol1/config/asevent_config.sh
source /myvol1/func/asevent_func.sh
tool=majiq


#make output directory 
mk_outdir $tool 
out=$wd/$output/${tool:-unspecific}-output #local variable for output folder


#handle SAM file
readsamfiles
for filename in $(cat $out/samlist)
do
	makebamfromsam $filename
done


#build custom config file, override if already present
config=$out/config.txt
echo "[info]" > $config
echo "readlen=$read_length" >> $config
#bam-directory depending on type of run
if [ $differential = 0 ]; then echo "bamdirs=$wd/$bamfolder" >> $config; else echo "bamdirs=$wd/$casefolder,$wd/$controlfolder" >> $config; fi
echo "genome=hg38" >> $config
echo "strandness=None" >> $config
echo "[experiments]" >> $config


#build config file depending on type of run
if [ $differential = 0 ]
then
	#list .bam files, remove last file extension (.bam), list them in one line with comma seperated, no comma after last file
	#(sed commands from here: https://unix.stackexchange.com/questions/313791/how-to-list-only-files-in-a-directory-separated-by-commas)
	bamlist=$(cd $wd/$bamfolder && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | xargs echo | sed 's/ /,/g')
	echo "ALLBAMS=$bamlist" >> $config

	echo building majiq reference ...
	majiq build $wd/$gff -c $config -j $ncores -o $out/build

	#get all .majiq files which were created with build 
	majiqlist=$(ls -1p $out/build/*.majiq | xargs echo)
	majiq psi $majiqlist -j $ncores -o $out/psi -n "ALLBAMS"
	cleaner
fi

if [ $differential = 1 ]
then
	# list .bam files from case & control folder, without file extension
	casefiles=$(cd $wd/$casefolder && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | xargs echo | sed 's/ /,/g')
	controlfiles=$(cd $wd/$controlfolder && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | xargs echo | sed 's/ /,/g')
	echo "CASE=$casefiles" >> $config
	echo "CONTROL=$controlfiles" >> $config

	majiq build $wd/$gff -c $config -j $ncores -o $out/build

	#list case/control files and add .majiq file extension as suffix and path of /build directory as prefix
	#(sed commands from here: https://askubuntu.com/questions/76808/how-do-i-use-variables-in-a-sed-command)
	caselist=$(cd $wd/$casefolder && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | sed "s|.*|$out/build/&.majiq|" | xargs echo)
	controllist=$(cd $wd/$controlfolder && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | sed "s|.*|$out/build/&.majiq|" | xargs echo)

	majiq deltapsi -grp1 $caselist -grp2 $controllist -n CASE CONTROL -j $ncores -o $out/deltapsi 
	cleaner
fi






