# test if fasta is available
test_fasta(){
	if ! test -f /$wd/$fasta; 
		then echo "File not found. Check the path for the ${fasta} fasta files: is it under ${wd}/${fasta}?";
		exit; fi
}
	
# test if gtf is available
test_gtf(){
	if ! test -f /$wd/$gtf; 
		then echo "File not found. Check the path for the ${gtf} gtf file: is it under ${wd}/${gtf}?";
		exit; fi
}

# make a list of fastq files to perform the mapping on and make it accessible
mk_fastqlist(){
find $fastqdir -name "*fastq" -nowarn -maxdepth 1| sed s/.fastq// | sed 's/.$//' | sort | uniq >/$wd/tmp/$tool-fastqlist
chmod 777 /$wd/tmp/$tool-fastqlist
}

#make the output directory and make it accessible
mk_outdir(){
mkdir -p /$wd/$out/
chmod -R 777 /$wd/$out/
}

#cleaning up: delete fastqlist from tmp folder
cleaner() {
# make output accessible
chmod -R 777 /$wd/$out/
rm /$wd/tmp/$tool-fastqlist
echo "script is done";
exit;
}