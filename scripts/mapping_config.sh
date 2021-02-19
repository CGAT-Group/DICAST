#####################
#     Parameters    #
#####################

### Uncomment to overwrite parameters set in the global config.sh ###
### which sets parameters for all tools (mapping and alternative splicing) ###

outname=$tool	# basename of output file (will usually be prefixed with the fastq file name and suffixed with .sam)


#################
#     Index     #
#################

recompute_index=false		#force index to be computed even if index with $indexname already exists
indexname=${fastaname}_index	#basename of index (without eg. .1.bt2 for bowtie index)
indexdir=$workdir/index/${tool:-unspecific}_index 	#directory of index

####################################
#     Tool specific parameters     #
####################################

# Some tools require chromosome-wise fasta-inputs:
contextmap_fastadir=$inputdir/fasta_chromosomes			#fasta !directory! for contextmap: chromosome wise fasta files
mapsplice_fastadir_mapping=$inputdir/fasta_chromosomes		#fasta !directory! for mapsplice mapping: chromosome wise fasta files
bowtie_fastadir=$inputdir/fasta_chromomes

