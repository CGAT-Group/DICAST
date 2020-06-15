# fasta = name of the fasta file to use.
# nthreads = number of threads to use
# inputdir = name of the input directory. Should contain all fasta and fastq files to use and only those to use. Type "" to use working directory ($wd) as input directory. 
#			Input files will be drawn from /$wd/$inputdir/

# tool = name of the tool. Defines output directory and file name
# wd = name of the mounted working directory. This directory should contain the input, output and index directories
# out = name of the output directory. The output will be saved to /$wd/$out/.
# index = name of the index file to use. A prebuild index has to be stored in /$wd/index/$tool-index/. 
#			If the index file is not found, the index will be build during runtime. This can take some time. 
#			To name the index after the corresponding fasta file use index = "$fasta.mmi"
#			Be careful to use the right index (build by the right tool and right fasta file)

### Only change if necessary:

# vars

tool=minimap2 
nthreads=64 

# directories

wd=myvol1
out=output/$tool-output 
inputdir=input
index=$fasta.mmi 

# filenames

#find .fa files in inputdir
fasta=$(find /$wd/$inputdir/ -maxdepth 1 -type f -name "*.fa" -printf "%f\n")
#gtf=Homo_sapiens.GRCh38.99.gtf
gtf=$(find /$wd/$inputdir/ -maxdepth 1 -type f -name "*.gtf" -printf "%f\n")
fastqdir=/$wd/$inputdir/