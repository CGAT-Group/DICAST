##########################
#     Test Parameters    #
##########################

# test if fasta is available
test_fasta(){
	if ! test -f $fasta;
		then echo "Fasta file not found. Is it under $fasta? \nCheck fasta variable config.sh and mapping_config.sh file! Parameters set in mapping_config.sh do overwrite config.sh!";
		exit;
	fi ;
}

# test if gtf is available (only tested if gtf is used for the mapping)
test_gtf(){
	if ! test -f $gtf;
		then echo "gtf file not found. Is it under $gtf? \nCheck gtf variable in config.sh and mapping_config.sh file! Parameters set in mapping_config.sh do overwrite config.sh!";
		exit;
	fi
}

####################
#     Fastqlist    #
####################

# make a list of fastq files to perform the mapping on and make it accessible
mk_fastqlist(){
	mkdir -p /tmp/
	find $fastqdir -name "*.fastq" -nowarn -maxdepth 2| sed s/.fastq$// | sed 's/.$//' | sort | uniq >/tmp/$tool-fastqlist
	while read -r line; do
		# mkdir -p $outdir/$(basename $(dirname $line))/
		mkdir -p $outdir/$(basename $(dirname $(dirname $line)))
		chmod -R 777 $outdir/
	done </tmp/$tool-fastqlist
}

######################
#     Directories    #
######################

#make the output directory and make it accessible
# mk_outdir(){
# 	mkdir -p $outdir/$(basename $(dirname $line))/
# 	chmod -R 777 $outdir/
# }

##################
#     Cleaner    #
##################

#cleaning up: delete fastqlist from tmp folder
cleaner() {
	# make output accessible
	chmod -R 777 $outdir/
	rm /tmp/$tool-fastqlist
	echo "Script is done.";
	exit;
}

#################
#    Logging    #
#################

echo_vars() {
	echo "############################################"
	echo "VARIABLE SETTINGS"
	echo -e """

tool:\t\t$tool

DIRECTORIES
workdir:\t$workdir
outdir: \t$outdir
inputdir:\t$inputdir
fastqdir:\t$fastqdir

FILES
fasta:\t\t$fasta
bowtie_fastadir:$bowtie_fastadir	# only for contextmap & mapsplice
gtfname:\t$gtfname	# only for hisat
gtf:\t\t$gtf # only for hisat mapsplice and star

INDEX
indexdir:\t$indexdir
indexname:\t$indexname
star_index:\t$star_index	# only for star

PARAMETERS
ncores: \t$ncores	# only for applicable tools
recompute_index:$recompute_index
"""
	echo "############################################"
	echo
}

start_logging() {
	mkdir -p $outdir/logs
	current_time=$(date "+%Y.%m.%d_%H:%M:%S")
	log_file=$outdir/logs/${tool}_${current_time}.log
	touch $log_file
	echo -e "\nlogs will be stored in ${log_file}\n"
	exec &> >(tee -a -i "$log_file")
	echo_vars
	#exec 2> >(tee -a -i "${log_file}")
	#exec >> "${log_file}"

}
