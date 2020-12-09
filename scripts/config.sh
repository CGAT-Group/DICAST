############################
#     Basic Parameters     #
############################

ncores=6                                                	#number of cores or threads the tool will use
workdir=/MOUNT/unify_test_MOUNT                         	#name of the base directory inside the Docker
outdir=$workdir/output/${tool:-unspecific}-output       	#name of the output directory; will be named after the specific tool that was used
read_length=76                                          	#length of reads inside fastq files

#############################
#     Input Directories     #
#############################


inputdir=$workdir/input         				#base directory for all needed input files
fastqdir=$inputdir/fastq        				#directory for fastqfiles
bamdir=$inputdir/bam            				#directory for bamfiles
samdir=$inputdir/sam            				#directory for samfiles
fastadir=$inputdir              				#directory for fastafile (might vary for specific tools -> see mapping or as-specific config file)
gtfdir=$inputdir                				#directory for gtffile
gffdir=$inputdir                				#directory for gfffile


############################
#     Input Parameters     #
############################

fastaname=Homo_sapiens.GRCh38.dna.primary_assembly.fa           #name of the genome reference file (fasta format), directory=$fastadir
gtfname=splicing_variants.gtf                                   #name of gtf reference file, directory=$gtffile
gffname=splicing_variants.gff3                                  #name of gff reference file, directory=$gfffile

fasta=${fastadir:-unspecific}/$fastaname                        #fasta full path
gtf=${gtfdir:-unspecific}/$gtfname                              #gtf full path
gff=${gffdir:-unspecific}/$gffname                              #gff full path

#################################
#     ASimulatoR Parameters     #
#################################

asimulator_inputdir=/nfs/home/students/tfaro/dockers/src/ASimulatoR/in
asimulator_outputdir=/nfs/home/students/tfaro/dockers/src/ASimulatoR/out




































