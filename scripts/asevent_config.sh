############################################
#  values for standard AS event detection  #
############################################

transcript=$inputdir/Transcripts.fa #$inputdir/custom_transcripts.fasta                  #fasta file for gene transcripts
	#gffread input/splicing_variants.gff3 -g input/Homo_sapiens.GRCh38.dna.primary_assembly.fa -w input/new_transcripts.fasta
star_alignment_files=$workdir/output/star-output                       #path to the folder containing star alignment files (*.SJ.)

#fastq files have to present as pairs; the next two fields contain the suffixes (including the file extension) all fastq-pairs present need to have (e.g. "_1.fastq" & "_2.fastq")
fastqpair1suffix="1.fastq"                                             #suffix for first file of fastq-pair
fastqpair2suffix="2.fastq"                                             #suffix for second file of fastq-pair
use_bam_input_files=1                                                   #determines what kind of input to use: 1 for bam files, 1 for fastq files

combine_events=1                                                        #For some tools, a unified output format is available. This sets wether or not the unification should combine multiple events into one event type.
