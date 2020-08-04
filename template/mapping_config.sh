# nthreads = number of threads to use

# wd = name of the mounted working directory. This directory should contain the input, output and index directories
# out = name of the output directory. The output will be saved to /$wd/$out/.
# inputdir = name of the input directory. Should contain all fasta and fastq files to use and only those to use. Type "" to use working directory ($wd) as input directory. 
#			Input files will be drawn from /$wd/$inputdir/
# index = basename of the index file to use. A prebuild index has to be stored in /$wd/index/$tool-index/. 
#			If the index file is not found, the index will be build during runtime. This can take some time. 
#			To name the index after the corresponding fasta file use index = "$fasta"
#			Be careful to use the right index (build by the right tool and right fasta file)
#			Note that some tools create multiple index files. In that case "$index" will be the basename of these files.
# fastqdir = name of the directory of the fastqfiles

# fasta = name of the fasta file to use. $(find /$wd/$inputdir/ -maxdepth 1 -type f -name "*.fa" -printf "%f\n") finds all .fa files in the input directory
# gtf = name of the gtf file to use. $(find /$wd/$inputdir/ -maxdepth 1 -type f -name "*.gtf" -printf "%f\n") finds all .gtf files in the input directory

### Only change if necessary:

# variables

nthreads=60
recompute_index=false

# directories

wd=myvol1
out=/$wd/output/${tool:-unspecific}-output
inputdir=/$wd/input
fastqdir=$inputdir/
indexdir=/$wd/index/${tool:-unspecific}-index

# filenames

#find .fa files in inputdir
fasta=$(find $inputdir -maxdepth 1 -type f -name "Homo*.fa" -printf "%f\n")
#gtf=Homo_sapiens.GRCh38.99.gtf
gtf=$(find $inputdir -maxdepth 1 -type f -name "*.gtf" -printf "%f\n")
index=$fasta
