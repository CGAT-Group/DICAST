#!/bin/bash


tool=majiq
source /MOUNT/scripts/config.sh
source /MOUNT/scripts/asevent_config.sh
source /MOUNT/scripts/asevent_func.sh

### logging ###
start_logging

#cleaning up
trap cleaner EXIT

#make output directory
mk_outdir


#handle SAM file
handlesamfiles $differential

#build custom config file, override if already present
config=$outdir/config.txt
echo "[info]" > $config
echo "readlen=$read_length" >> $config
#bam-directory depending on type of run
echo "bamdirs=$casebam,$controlbam" >> $config
echo "genome=hg38" >> $config
echo "strandness=None" >> $config
echo "[experiments]" >> $config


#build config file depending on type of run
if [ $differential = 0 ]
then
	#list .bam files, remove last file extension (.bam), list them in one line with comma seperated, no comma after last file
	#(sed commands from here: https://unix.stackexchange.com/questions/313791/how-to-list-only-files-in-a-directory-separated-by-commas)
	#bamlist=$(cd $controlbam && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | xargs echo | sed 's/ /,/g')
	for filename in $(cat /tmp/controlbamlist)
	do
		# build sperate config file for each BAM file
		basename=$(basename $filename)
		majiq_basename=$(basename -s .bam $filename)
		outdir_name=$(basename -s .bam $filename)_output
		mkdir -p $outdir/$outdir_name
		config=$outdir/$outdir_name/config.txt
		echo "[info]" > $config
		echo "readlen=$read_length" >> $config
		#bam-directory depending on type of run
		echo "bamdirs=$controlbam" >> $config
		echo "genome=hg38" >> $config
		echo "strandness=None" >> $config
		echo "[experiments]" >> $config
		echo "BAM=$majiq_basename" >> $config

		echo "building MAJIQ reference ..."
		majiq build $gff -c $config -j $ncores -o $outdir/$outdir_name/build
		wait

		#get all .majiq files which were created with build
	        majiqlist=$(ls -1p $outdir/$outdir_name/build/*.majiq | xargs echo)

		majiq psi $majiqlist -j $ncores -o $outdir/$outdir_name/psi -n "BAM"

		# create voila.tsv outputfiles
		voila tsv $outdir/$outdir_name/build/splicegraph.sql $outdir/$outdir_name/psi/*.voila -f $outdir/$outdir_name/voila.tsv
		wait

		echo "Running $tool unificiation..."

		echo "Looking for whippet files in $outdir/$outdir_name"
		unified_outdir_name="${outdir}/${outdir_name}_${tool}_dicast_unify"
		mkdir -p $unified_outdir_name
		echo "Saving unified output to $unified_outdir_name"

		if [[ -f "$outdir/$outdir_name/psi/BAM.psi.tsv" && -f "$outdir/$outdir_name/voila.tsv" ]]; then

			# Save all majiq output files to the same new tmp directory
			uni_tmp="/tmp/unification_tmpdir"
			mkdir -p $uni_tmp
			ln -sf $outdir/$outdir_name/psi/BAM.psi.tsv $uni_tmp/BAM.psi.tsv
			ln -sf $outdir/$outdir_name/voila.tsv $uni_tmp/voila.tsv
			mkdir -p $unified_outdir_name

			anno_file="$workdir/src/ASimulatoR/out/event_annotation.tsv"
			stats_file="${unified_outdir_name}/${outdir_name}_${tool}_dicast_unify_comparison.txt"

			if [ $combine_events = 0 ];
			then
				python3 /MOUNT/scripts/unified_output/output_transformer.py create -m $uni_tmp -out $unified_outdir_name -gtf $gtf

				if [[ -f "$anno_file" ]];
				then
					echo "Running unified comparison..."
					python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_${tool}_dicast_unify.out -gtf $gtf -stats $stats_file -s -t 0
				fi
			else
				python3 /MOUNT/scripts/unified_output/output_transformer.py create -m $uni_tmp -out $unified_outdir_name -gtf $gtf -comb
				if [[ -f "$anno_file" ]];
				then
					echo "Running unified comparison..."
					python3 /MOUNT/scripts/unified_output/output_transformer.py compare -a $anno_file -c ${unified_outdir_name}/${outdir_name}_${tool}_dicast_unify.out -gtf $gtf -stats $stats_file -s -t 0 -comb
				fi
			fi
		else
			echo "Couldn't find necessary input files for unification: $outdir/$outdir_name/psi/BAM.psi.tsv and $outdir/$outdir_name/voila.tsv"
		fi

		echo "Finished $tool unification for ${outdir_name}."
	done
	cleaner
fi

if [ $differential = 1 ]
then
	#build custom config file, override if already present
	config=$outdir/config.txt
	echo "[info]" > $config
	echo "readlen=$read_length" >> $config
	#bam-directory depending on type of run
	echo "bamdirs=$casebam,$controlbam" >> $config
	echo "genome=hg38" >> $config
	echo "strandness=None" >> $config
	echo "[experiments]" >> $config


	# list .bam files from case & control folder, without file extension
	casefiles=$(cd $casebam && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | xargs echo | sed 's/ /,/g')
	controlfiles=$(cd $controlbam && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | xargs echo | sed 's/ /,/g')
	echo "CASE=$casefiles" >> $config
	echo "CONTROL=$controlfiles" >> $config

	majiq build $gff -c $config -j $ncores -o $outdir/build

	#list case/control files and add .majiq file extension as suffix and path of /build directory as prefix
	#(sed commands from here: https://askubuntu.com/questions/76808/how-do-i-use-variables-in-a-sed-command)
	caselist=$(cd $casebam && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | sed "s|.*|$outdir/build/&.majiq|" | xargs echo)
	controllist=$(cd $controlbam && ls -1p *.bam | grep -Po '.*(?=\.)' | grep -v / | sed "s|.*|$outdir/build/&.majiq|" | xargs echo)

	majiq deltapsi -grp1 $caselist -grp2 $controllist -n CASE CONTROL -j $ncores -o $outdir/deltapsi

	#create voila.tsv outputfile
	voila tsv $outdir/build/splicegraph.sql $outdir/deltapsi/*.voila -f $outdir/voila.tsv
	cleaner
fi
