############################################
#  values for standard AS event detection  #
############################################

transcript=input_asevent/custom_transcripts.fasta			#fasta file for gene transcripts

#fastq files have to present as pairs; the next two fields contain the suffixes (including the file extension) all fastq-pairs present need to have (e.g. "_1.fastq" & "_2.fastq")
fastqpair1suffix="_1.fastq"						#suffix for first file of fastq-pair 
fastqpair2suffix="_2.fastq"						#suffix for second file of fastq-pair


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
