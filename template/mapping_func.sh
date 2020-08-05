# test if fasta is available
test_fasta(){
	if ! test -f $inputdir/$fasta; 
		then echo "File not found. Check the path for the ${fasta} fasta files: is it under $inputdir/$fasta?";
		exit; fi
}
	
# test if gtf is available (only tested if gtf is used for the mapping)
test_gtf(){
	if ! test -f $inputdir/$gtf; 
		then echo "File not found. Check the path for the ${gtf} gtf file: is it under $inputdir/$gtf?";
		exit; fi
}

# make a list of fastq files to perform the mapping on and make it accessible
mk_fastqlist(){
	find $fastqdir -name "*fastq" -nowarn -maxdepth 1| sed s/.fastq// | sed 's/.$//' | sort | uniq > /tmp/$tool-fastqlist
	chmod 777 $out/$tool-fastqlist
}

#make the output directory and make it accessible
mk_outdir(){
	mkdir -p $out/
	chmod -R 777 $out/
}

#cleaning up: delete fastqlist from tmp folder
cleaner() {
	# make output accessible
	chmod -R 777 $out/
	rm $out/$tool-fastqlist
	echo "Script is done.";
	exit;
}
