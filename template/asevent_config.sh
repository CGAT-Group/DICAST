###################
#  basic values   #
###################

wd=/myvol1                                                              #name of mounted folder in docker container (don't change, unless you know what you are doing); all other path variables will start from here
output=output                                                           #name of output folder (the tool name will be added to this name with "-$tool")
ncores=60                                                               #number of cores to use


############################################
#  values for standard AS event detection  #
############################################

gtf=input_asevent/splicing_variants.gtf					#name of gtf-file; if not in $wd, give full path here
#gtf=input_asevent/transcripts.gtf
gff=input_asevent/splicing_variants.gff3				#name of gff3 file of reference
bamfolder=input_asevent/bamfolder                                       #path of BAM-folder; always expecting BAM+BAI file
fasta=input_asevent/genome.fa						#genome reference file (used by IRFinder, KisSplice and ASGAL)
transcript=input_asevent/custom_transcripts.fasta			#fasta file for gene transcripts
fastqfolder=input_asevent/fastqfolder					#path of fastq-folder; all files which want to be used have to end with .fastq

#fastq files have to present as pairs; the next two fields contain the suffixes (including the file extension) all fastq-pairs present need to have (e.g. "_1.fastq" & "_2.fastq")
fastqpair1suffix="_1.fastq"						#suffix for first file of fastq-pair 
fastqpair2suffix="_2.fastq"						#suffix for second file of fastq-pair

read_length=76								#length of mapped reads 
use_bam_input_files=1							#determines what kind of input to use: 1 for bam files, 0 for fastq files

#####################################
# values for differential analysis  #
#####################################

differential=1 								#1: tools which can calculate differential splicing, will use it; 0: only AS event detection
casefolder=input_asevent/casefolder					#path of folder with BAMs used as case for DS analysis; needs to be filled when differential=1; filled like bamfolder
casefastq=input_asevent/casefastq					#path of folder with fastq files used as case for DS analysis
caseprefix="subsample_0"						#all files in the casefastq folder must have this by the user specified prefix
controlfolder=input_asevent/controlfolder				#path of folder with BAMs used as control for DS analysis; needs to be filled when differential=1; filled like bamfolder
controlfastq=input_asevent/controlfastq					#path og folder with fastq files used as control for DS analysis
controlprefix="subsample_1"						#add files in the controlfastq folder must have this by the user specified prefix
