<<<<<<< HEAD
#!/bin/bash

### Tool-specific variables ###
tool=segemehl

# use confic and function file
source /myvol1/config/mapping_config.sh
source /myvol1/func/mapping_func.sh


### Tool-specific functions ###

# Unpaired mapping command: second attempt, used if paired end mapping failes
second_attempt() {
	# tag outputs with this flag to name it per fastqfile 	"${line##*/}"
	echo "paired mapping failed for ${line}. Try unpaired mapping."
	#Parameters
	# -d 	--database [<file>] 	list of filename(s) of fasta database sequence(s)
	# -q 	--query 	filename of query sequences
	# -p 	--mate 	filename of paired end sequences
	# -i 	--index 	filename of database index
	#X -j 	--index2 	filename of second database index
	#X -x 	--generate 	generate db index and store to disk
	#X -y 	--generate2 	generate second db index and store to disk
	#X -G 	--readgroupfile 	filename to read @RG header (default:none)
	#X -g 	--readgroupid 	read group id (default:none)
	# -t 	--threads 	start threads (default:1)
	#X -f 	--fullname 	write full fastq/a name to SAM
	segemehl.x \
		-d $inputdir/$fasta \
		-q ${line}?.fastq \
		-i $indexdir/$index \
		-t $nthreads \
		-o $out/${line##*/}${tool}.sam
}

# make index directory and build index if index was not found
build_index() {
	mkdir -p $indexdir
	echo "compute index"
	segemehl.x -x $indexdir/$index -d $(ls $inputdir/$fasta)
	chmod -R 777 $indexdir
	echo "Index is now saved under /$wd/index/$tool-index/$index.idx"
}

### START here ############################################################################

# test filepaths
test_fasta

# Build Genome index if not already available
if $recompute_index; then build_index; else if ! test -f $indexdir/$index; then build_index; fi fi

#make list of fastq files
mk_fastqlist

#make output directories
mk_outdir

### Start mapping ###

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do
	#First attempt: Paired end mapping
	#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
	#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")
	
	# Parameters (X are not used)
	# -d 	--database [<file>] 	list of filename(s) of fasta database sequence(s)
	# -q 	--query 	filename of query sequences
	# -p 	--mate 	filename of paired end sequences
	# -i 	--index 	filename of database index
	#X -j 	--index2 	filename of second database index
	#X -x 	--generate 	generate db index and store to disk
	#X -y 	--generate2 	generate second db index and store to disk
	#X -G 	--readgroupfile 	filename to read @RG header (default:none)
	#X -g 	--readgroupid 	read group id (default:none)
	# -t 	--threads 	start threads (default:1)
	#X -f 	--fullname 	write full fastq/a name to SAM

segemehl.x \
	-d $inputdir/$fasta \
	-q ${line}1.fastq  \
	-p ${line}2.fastq \
	-i $indexdir/$index \
	-t $nthreads \
	-o $out/${line##*/}${tool}.sam

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </tmp/$tool-fastqlist

# wait for all processes to end
wait
cleaner
||||||| merged common ancestors
=======
#!/bin/bash

### Needed Functions

#Unpaired mapping command
#...tag outputs with this flag to name it per fastqfile 	"${line##*/}"
#...address for all gtf files are 				$(find /myvol1/ -name "*.gtf")
second_attempt() {
echo "paired mapping failed for ${line}. Try unpaired mapping."
segemehl.x \
	-d /myvol1/${fasta} \
	-q ${line}.fastq \
	-i /myvol1/index/${tool}-index/${fasta}.idx \
	-t ${nthreads} \
	-o /myvol1/output/"$tool"-output/${line##*/}${tool}.sam
}

build_index() {
#Build Genome index
# find gtf and da files with these  $(find /myvol1/ -name "*.fa")  $(find /myvol1/ -name "*.gtf")
mkdir -p /myvol1/index/${tool}-index
echo "compute index"
segemehl.x -x /myvol1/index/${tool}-index/${fasta}.idx -d $(ls /myvol1/${fasta})
chmod -R 777 /myvol1/index/${tool}-index/
}

#cleaning up
cleaner() {
rm /myvol1/tmp/"$tool"-fastqlist
echo "script is done";exit;
}

#START here: Make a list of fastq files ############################################################################
tool=segemehl
nthreads=64
fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa

#test filepaths for fasta and indexing
if ! test -f "/myvol1/${fasta}"; 
	then echo "check the path for the ${fasta} fasta files: is it under <mounted folder>/${fasta}?";
	cleaner;
	exit; fi
#if ! test -f "/myvol1/index/${tool}-index/${fasta}.idx"; then echo "would build index now"; fi
if ! test -f "/myvol1/index/${tool}-index/${fasta}.idx"; then build_index; fi

find /myvol1/ -name "*fastq" -nowarn -maxdepth 1| sed s/.fastq// | sed 's/.$//' | sort | uniq >/myvol1/tmp/"$tool"-fastqlist
chmod 777 /myvol1/tmp/"$tool"-fastqlist

#make output directories
mkdir -p /myvol1/output/"$tool"-output/
chmod -R 777 /myvol1/output/"$tool"-output/

echo "compute ${tool} mapping..."
#Iterate list with paired end map command first
while read -r line; do

#First attempt: Paired end mapping
#...tag outputs with this flag to name it per fastqfile         "${line##*/}"
#...address for all gtf files are                               $(find /myvol1/ -name "*.gtf")

# -d 	--database [<file>] 	list of filename(s) of fasta database sequence(s)
# -q 	--query 	filename of query sequences
# -p 	--mate 	filename of paired end sequences
# -i 	--index 	filename of database index
#X -j 	--index2 	filename of second database index
#X -x 	--generate 	generate db index and store to disk
#X -y 	--generate2 	generate second db index and store to disk
#X -G 	--readgroupfile 	filename to read @RG header (default:none)
#X -g 	--readgroupid 	read group id (default:none)
# -t 	--threads 	start threads (default:1)
#X -f 	--fullname 	write full fastq/a name to SAM

segemehl.x \
	-d /myvol1/${fasta} \
	-q ${line}1.fastq  \
	-p ${line}2.fastq \
	-i /myvol1/index/${tool}-index/${fasta}.idx \
	-t ${nthreads} \
	-o /myvol1/output/"$tool"-output/${line##*/}${tool}.sam

#If paired end mapping fails, run unpaired mapping.
trap 'second_attempt $line' ERR
done </myvol1/tmp/"$tool"-fastqlist

#wait for all processes to end
wait
#make output accessible
chmod -R 777 /myvol1/output/"$tool"-output/
cleaner
>>>>>>> origin/minimap_entry
