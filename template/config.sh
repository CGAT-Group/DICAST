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
fastqdir=$workdir/$inputdir/fastq					#directory for fastqfiles
bamdir=$workdir/$inputdir/bamdir					#directory for bamfiles
samdir=$workdir/$inputdir/samdir					#directory for samfiles
fasta=$workdir/$inputdir/fasta						#name of the genome reference file (fasta format)
gtf=$workdir/$inputdir/gtf						#name of gtf reference file
gff=$workdir/$inputdir/gff						#name of gff reference file
indexdir=$workdir/$inputdir/index					#directory for STAR index
indexname=$fasta							#name of the STAR index



