############################################
#  values for standard AS event detection  #
############################################

transcript=$fasta #$inputdir/custom_transcripts.fasta                  #fasta file for gene transcripts
star_alignment_files=$workdir/output/star-output                       #path to the folder containing star alignment files (*.SJ.)

#fastq files have to present as pairs; the next two fields contain the suffixes (including the file extension) all fastq-pairs present need to have (e.g. "_1.fastq" & "_2.fastq")
fastqpair1suffix="_1.fastq"                                             #suffix for first file of fastq-pair
fastqpair2suffix="_2.fastq"                                             #suffix for second file of fastq-pair
use_bam_input_files=0                                                   #determines what kind of input to use: 1 for bam files, 0 for fastq files

#####################################
# values for differential analysis  #
#####################################

differential=0                                                          #1: tools which can calculate differential splicing, will use it; 0: only AS event detection

casebam=$casefolder/bamdir                                              #path of folder with BAMs used as case for DS analysis; needs to be filled when differential=1; filled like bamfolder
casefastq=$casefolder/fastqdir                                          #path of folder with fastq files used as case for DS analysis
caseprefix="sample_01"                                                #all files in the casefastq folder must have this by the user specified prefix

controlbam=$controlfolder/bamdir                                        #path of folder with BAMs used as control for DS analysis; needs to be filled when differential=1; filled like bamfolder
controlfastq=$controlfolder/fastqdir                                    #path of folder with fastq files used as control for DS analysis
controlprefix="sample_01"                                             #add files in the controlfastq folder must have this by the user specified prefix

#########################################
# Paths to renamed files for PSI-Sigma  #
#########################################
renamed_casebam=$casebam/renamed_psi_sigma/                             #like casebam, but file names need special formatting, see documentation
renamed_controlbam=$controlbam/renamed_psi_sigma/                       #like controlbam, but file names need special formatting, see documentation
renamed_star_alignment_files=$star_alignment_files/renamed_psi_sigma/   #like star_alignment_files, but file names need special formatting, see documentation
