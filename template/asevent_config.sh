#"tool" variable (tool-name) will be set in each ENTRYPOINT script seperatly

wd=/myvol1								#name of mounted folder in docker container (don't change, unless you know what you are doing); all other path variables will start from here
output=output								#name of output folder (the tool name will be added to this name with "-$tool")
gtf=splicing_variants.gtf						#name of gtf-file; if not in $wd, give full path here
bamfile=Star_mapped_subsample_01_100K_Aligned_sorted_filtered.out.bam	#name of BAM file; has to be within bamfolder
fastqfolder=input_asevent/fastq						#path of fastq-folder
bamfolder=subsampled/star-output					#path of BAM-folder
read_length=76								#length of mapped reads in BAM file
ncores=60								#number of cores to use
differential=0 								#1: tools which can calculate differential splicing, will use it; 0: only AS event detection
casefolder=								#path of folder with BAMs used as case for DS analysis; needs to be filled when differential=1
controlfolder= 								#path of folder with BAMs used as control for DS analysis; needs to be filled when differential=1

