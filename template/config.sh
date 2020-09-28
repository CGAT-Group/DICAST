############################
#     Basic Parameters     #
############################

ncores=6										#number of cores the tool will use
workdir=/MOUNT								#name of the base directory inside the Docker
outdir=$workdir/output/${tool:-unspecific}-output			#name of the output directory; will be named after the specific tool that was used
read_length=76								#length of reads inside fastq files

############################
#     Input Parameters     #
############################

fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa 		#name of the genome reference file (fasta format), directory=$fastadir
gtf=splicing_variants.gtf										#name of gtf reference file, directory=$gtffile
gff=splicing_variants.gff3									#name of gff reference file, directory=$gfffile


#############################
#     Input Directories     #
#############################


inputdir=$workdir/input		#base directory for all needed input files
fastqdir=$inputdir/fastq		#directory for fastqfiles
bamdir=$inputdir/bam			#directory for bamfiles
samdir=$inputdir/sam			#directory for samfiles
fastadir=$inputdir				#directory for fastafile (might vary for specific tools -> see mapping or as-specific config file)
gtfdir=$inputdir				#directory for gtffile
gffdir=$inputdir				#directory for gfffile




