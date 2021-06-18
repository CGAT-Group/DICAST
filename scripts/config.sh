############################
#     Basic Parameters     #
############################

ncores=16                                                	#number of cores or threads the tool will use
workdir=/MOUNT			                         	#name of the base directory inside the Docker
outdir=$workdir/output/${tool:-unspecific}-output       	#name of the output directory; will be named after the specific tool that was used
read_length=76                                          	#length of reads inside fastq files

#############################
#     Input Directories     #
#############################


inputdir=$workdir/input
controlfolder=$inputdir/controldir         			#base directory for all needed input files (when no differential comparison, control inputs when differential AS Event Detection)
casefolder=$inputdir/casedir					#base directory for only case files (for AS Event detection)
fastqdir=$controlfolder/fastqdir       				#directory for fastqfiles
bamdir=$controlfolder/bamdir           				#directory for bamfiles
samdir=$controlfolder/bamdir           				#directory for samfiles
fastadir=$inputdir              				#directory for fastafile (might vary for specific tools -> see mapping or as-specific config file)
gtfdir=$inputdir                				#directory for gtffile
gffdir=$inputdir                				#directory for gfffile
star_index=$workdir/index/star_index                                   #folder containing a star index built with the $gtf and $fasta files (used by: IRFinder, KisSplice, rMATS)


############################
#     Input Parameters     #
############################

asimulator_gtf=Homo_sapiens.GRCh38.104.gtf			#name of the GTF file used to generate simulated data within ASimulatoR R library.
fastaname=Homo_sapiens.GRCh38.dna.primary_assembly.fa           #name of the genome reference file (fasta format), directory=$fastadir
gtfname=ASimulatoR.gtf        	       	                #name of gtf reference file, directory=$gtffile; set to ASimulatoR_gtf.gtf, when ASimulator is true
gffname=ASimulatoR.gff3					#set to ASimulatoR_gff.gff3, when ASimulator is true


fasta=${fastadir}/$fastaname                        #fasta full path
gtf=${gtfdir}/$gtfname                              #gtf full path
gff=${gffdir}/$gffname                              #gff full path