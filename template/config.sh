############################
#     Basic Parameters     #
############################

ncores=6								#number of cores the tool will use
workdir=/MOUNT								#name of the base directory inside the Docker
outdir=$workdir/output/${tool:-unspecific}-output			#name of the output directory; will be named after the specific tool that was used
read_length=76								#length of reads inside fastq files


############################
#     Input Parameters     #
############################


inputdir=$workdir/input							#base directory for all needed input files
fastqdir=$inputdir/fastq						#directory for fastqfiles
bamdir=$inputdir/bamdir							#directory for bamfiles
samdir=$inputdir/samdir							#directory for samfiles
fasta=$inputdir/fasta							#name of the genome reference file (fasta format)
gtf=$inputdir/gtf							#name of gtf reference file
gff=$inputdir/gff							#name of gff reference file
indexdir=$inputdir/index						#directory for STAR index
indexname=$fasta							#name of the STAR index



