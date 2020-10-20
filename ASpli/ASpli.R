<<<<<<< HEAD
suppressMessages(library(optparse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(ASpli))

#Parsing parameters
option_list = list(make_option(c("--gtf"), type='character', default = NULL, help="the annotation gtf file", metavar='character'),
                   make_option(c("--cores"), type='integer', default = 1, help="the number of cores", metavar='character'),
                   make_option(c("--readLength"), type='integer', default = NULL, help="the read length", metavar='character'),
		   make_option(c("--out"), type ='character',default = NULL, help="output directory",metavar='character'),
		   make_option(c("--differential"),type='integer',help="1 for differential analysis, needs at least two BAMs in config (case&control); 0 for AS event detection",metavar='character'),
		   make_option(c("--workdir"),type='character',help="the working drectory (all of the used folders have to be inside here)",metavar="character"),
		   make_option(c("--bamfolder"),type='character',default=NULL,help="folder containing BAM file used for AS event detection; only use this paramter if --differential parameter is NOT used",metavar='character'),
		   make_option(c("--bamfile"),type='character',default=NULL,help="Name of BAM file to use for AS event detection"),
		   make_option(c("--casefolder"),type='character',default=NULL,help="folder containing all BAM-files, which will be labeled CASE for differential splicing analysis; --differential parameter needs to be used", metavar='character'),
		   make_option(c("--controlfolder"),type='character',default=NULL,help="folder containing all BAM-files, which will be labeled CONTROL for differential splicing analysis; --differential parameter needs to be used", metavar='character'))



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

differential <- ifelse(opt$differential,1,0) #if differential==1 --> TRUE

if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide an annotation gtf file", call.=FALSE)
}
if (is.null(opt$readLength)){
  print_help(opt_parser)
  stop("Please provide reads length", call.=FALSE)
}
if(is.null(opt$out)){
  print_help(opt_parser)
  stop("Please provide an output directory", call.=FALSE)
}
if(is.null(opt$workdir)){
  print_help(opt_parser)
  stop("Please provide the working directory", call.=FALSE)
}
if(differential){
  if(is.null(opt$casefolder)) stop("Please provide a case-folder")
  if(is.null(opt$controlfolder)) stop ("Please provide control-folder")
}else{
  if (is.null(opt$bamfolder)) stop ("Please provide bamfolder")
  #if (is.null(opt$bamfile)) stop ("Please provide bamfile")
}


cores=opt$cores
readLength=opt$readLength
### build targtes-dataframe ###

if(differential){
  #get files from case & control folders with full path, but only BAM files, not the index!
  cases <- list.files(paste0(opt$workdir,"/",opt$casefolder),full.names=T, pattern="\\.bam$")
  controls <- list.files(paste0(opt$workdir,"/",opt$controlfolder),full.names=T,pattern="\\.bam$")

  targets <- data.frame(row.names=paste0("Sample_",c(1:(length(cases)+length(controls)))),bam=c(cases,controls),genotype=c(rep("Case",length(cases)),rep("Control",length(controls))))
  #targets$bam <- paste0(opt$workdir,targets$bam) #add workdir to each filepath of BAMfiles
}else{
  asevent_bams <- list.files(paste0(opt$workdir,"/",opt$bamfolder),full.names=T,pattern="\\.bam$")
  #asevent_bams <- basename(opt$bamfile)
  targets <- data.frame(row.names=paste0("Sample_",c(1:length(asevent_bams))),bam=asevent_bams,genotype=paste0("BAM",c(1:length(asevent_bams))))
}
print(paste0("Building TxDB with ",opt$gtf,"..."))
#Building a TxDb
TxDb <- makeTxDbFromGFF(file=opt$gtf, format="gtf") #the gtf file

# extract features from annotation
features <- binGenome(TxDb,logTo=paste0(opt$out,"/ASpli_binFeatures.log"))

print("loading BAM-file(s)...")
#loading bam
bam <- loadBAM(targets,cores=cores)

print("counting reads...")
#Reads counting
counts <- readCounts (features, bam, targets, cores = cores, readLength = readLength, maxISize = 50000, minAnchor = 10 ) #cores, readLength
writeCounts(counts=counts, output.dir = paste0(opt$out,"/ASpli_counts"))

print("discovering as events...")
#AS discovery

as <- AsDiscover( counts, targets, features, bam, readLength=readLength, threshold = 5, cores = cores)
writeAS(as=as, output.dir=paste0(opt$out,"/ASpli_as"))

#Differential Analysis
if(differential){
	print("differential analysis starting...")
	du <- DUreport( counts, targets)
	writeDU(du, output.dir=paste0(opt$out,"/ASpli_du"))
}
||||||| merged common ancestors
=======
suppressMessages(library(optparse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(ASpli))

#Parsing parameters
option_list = list(make_option(c("--gtf"), type='character', default = NULL, help="the annotation gtf file", metavar='character'),
                   make_option(c("--cores"), type='integer', default = 1, help="the number of cores", metavar='character'),
                   make_option(c("--readLength"), type='integer', default = NULL, help="the read length", metavar='character'),
		   make_option(c("--out"), type ='character',default = NULL, help="output directory",metavar='character'),
		   make_option(c("--differential"),type='integer',help="1 for differential analysis, needs at least two BAMs in config (case&control); 0 for AS event detection",metavar='character'), 
		   make_option(c("--bamfolder"),type='character',default=NULL,help="folder containing BAM file used for AS event detection; only use this paramter if --differential parameter is NOT used",metavar='character'),
		   make_option(c("--bamfile"),type='character',default=NULL,help="Name of BAM file to use for AS event detection"),
		   make_option(c("--casefolder"),type='character',default=NULL,help="folder containing all BAM-files, which will be labeled CASE for differential splicing analysis; --differential parameter needs to be used", metavar='character'),
		   make_option(c("--controlfolder"),type='character',default=NULL,help="folder containing all BAM-files, which will be labeled CONTROL for differential splicing analysis; --differential parameter needs to be used", metavar='character'))



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

differential <- ifelse(opt$differential,1,0) #if differential==1 --> TRUE

if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide an annotation gtf file", call.=FALSE)
}
if (is.null(opt$readLength)){
  print_help(opt_parser)
  stop("Please provide reads length", call.=FALSE)
}
if(is.null(opt$out)){
  print_help(opt_parser)
  stop("Please provide an output directory", call.=FALSE)
}
if(differential){
  if(is.null(opt$casefolder)) stop("Please provide a case-folder")
  if(is.null(opt$controlfolder)) stop ("Please provide control-folder")
}else{
  if (is.null(opt$bamfolder)) stop ("Please provide bamfolder")
  #if (is.null(opt$bamfile)) stop ("Please provide bamfile")
}


cores=opt$cores
readLength=opt$readLength
### build targtes-dataframe ###

if(differential){
  #get files from case & control folders with full path, but only BAM files, not the index!
  cases <- list.files(opt$casefolder,full.names=T, pattern="\\.bam$")
  controls <- list.files(opt$controlfolder,full.names=T,pattern="\\.bam$")
	
  targets <- data.frame(row.names=paste0("Sample_",c(1:(length(cases)+length(controls)))),bam=c(cases,controls),genotype=c(rep("Case",length(cases)),rep("Control",length(controls))))
  #targets$bam <- paste0(opt$workdir,targets$bam) #add workdir to each filepath of BAMfiles
}else{
  asevent_bams <- list.files(opt$bamfolder,full.names=T,pattern="\\.bam$")
  #asevent_bams <- basename(opt$bamfile)
  targets <- data.frame(row.names=paste0("Sample_",c(1:length(asevent_bams))),bam=asevent_bams,genotype=paste0("BAM",c(1:length(asevent_bams))))
}
print(paste0("Building TxDB with ",opt$gtf,"..."))
#Building a TxDb
TxDb <- makeTxDbFromGFF(file=opt$gtf, format="gtf") #the gtf file

# extract features from annotation
features <- binGenome(TxDb,logTo=paste0(opt$out,"/ASpli_binFeatures.log"))

print("loading BAM-file(s)...")
#loading bam
bam <- loadBAM(targets,cores=cores)

print("counting reads...")
#Reads counting
counts <- readCounts (features, bam, targets, cores = cores, readLength = readLength, maxISize = 50000, minAnchor = 10 ) #cores, readLength
writeCounts(counts=counts, output.dir = paste0(opt$out,"/ASpli_counts"))

print("discovering as events...")
#AS discovery

as <- AsDiscover( counts, targets, features, bam, readLength=readLength, threshold = 5, cores = cores)
writeAS(as=as, output.dir=paste0(opt$out,"/ASpli_as"))

#Differential Analysis
if(differential){
	print("differential analysis starting...")
	du <- DUreport( counts, targets)
	writeDU(du, output.dir=paste0(opt$out,"/ASpli_du"))
}
>>>>>>> unify-config
