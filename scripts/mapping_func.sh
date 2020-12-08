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
}

######################
#     Directories    #
######################

#make the output directory and make it accessible
mk_outdir(){
	mkdir -p $outdir/{case,control}
	chmod -R 777 $outdir/
}

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
