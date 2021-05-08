suppressMessages(library(EventPointer))
suppressMessages(library(optparse))

#Parsing arguments
option_list = list(make_option(c("--gtf"), type='character', default = NULL, help="the annotation gtf file", metavar='character'),
                   make_option(c("--bamfolder"), type='character', default = NULL, help="the path to bam file-folder", metavar='character'),
		   make_option(c("--bamfile"), type='character',default=NULL,help="Name of BAM-file for AS event detection"),
		   make_option(c("--differential"),type='integer',default=0,help="0:only AS event detection using file from --bamfile;1: differential analysis using files from --casefolder & --controlfolder"),
		   make_option(c("--casefolder"),type='character',default=NULL,help="path to folder with BAMs considered case (only use with --differential 1)"),
		   make_option(c("--controlfolder"),type='character',default=NULL,help="path to folder with BAMs considered control (only use with --differential 1)"),
		   make_option(c("--combined"),type="character",default=NULL,help="path to folder where files from case & control are located"),
		   make_option(c("--output"),type='character',default=NULL, help="output folder"),
		   make_option(c("--cores"), type='integer',default=1,help="Number of cores",metavar="integer"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

differential <- ifelse(opt$differential,1,0)  # TRUE if --differential == 1

if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide an annotation gtf file", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop ("Please provide an output directory")
}
if(differential){
  if(is.null(opt$casefolder)) stop("Please provide a case-folder")
  if(is.null(opt$controlfolder)) stop ("Please provide a control-folder")
  if(is.null(opt$combined)) stop("Please provide combined folder")
}else{
  if(is.null(opt$bamfolder)) stop("Please provide a bamfolder")
}


cores = as.numeric(opt$cores)
#if(getwd() == opt$workdir){opt$workdir <- ""}  #set opt$workdir empty if it is the current directory

PathToGTF<- opt$gtf

if(differential){

  # read combined folder
  Samples <- basename(list.files(opt$combined,pattern="\\.bam$"))
  PathToSamples <- opt$combined

  print("Preparing BAM files ...")

  SG_RNASeq<-PrepareBam_EP(Samples=bam_file,
                           SamplePath=PathToSamples,
                           Ref_Transc="GTF",
                           fileTransc=PathToGTF,
                           cores=cores)

  print("Looking for AS events...")
  TxtPath<-paste0(opt$output)
  AllEvents_RNASeq<-EventDetection(SG_RNASeq, cores=cores, Path=TxtPath)

  print("Starting differential analysis...")
  #number_of_one = opt$one
  #number_of_two = opt$two
  number_of_one <- length(list.files(opt$case))/2  #number of files in case folder; /2 to account for .bai files
  number_of_two <- length(list.files(opt$control))/2  #number of files in control folder; /2 to account for .bai file
  Cmatrix<-t(t(c(0,1))) #for two conditions
  Dmatrix <- matrix(c(rep(1, length(Samples)), rep(0, number_of_one), rep(1, number_of_two)), ncol=2,byrow=FALSE)
  print(Dmatrix)
  Events <- EventPointer_RNASeq(AllEvents_RNASeq,Dmatrix,Cmatrix,Statistic="LogFC",PSI=TRUE)
  write.table(Events, paste0(opt$output,"/EP_DAS.txt"), sep='\t')

}else{

  samplelist <- basename(list.files(opt$bamfolder,pattern="\\.bam$"))
  PathToSamples <- opt$bamfolder
  print(paste("Run Eventpointer for found samples:", samplelist))

  for (i in range(1:length(samplelist))){
      print(cat("Samples iteration", i, ":", samplelist, samplelist[i], "\n"))
      tryCatch({
      bam_file = samplelist[i]
      print(paste("Preparing BAM file",i,":",bam_file,"..." ))

      TxtPath<-paste0(opt$output,"/",bam_file,"_output")
      if(!dir.exists(TxtPath)){
        dir.create(TxtPath)
      }

      SG_RNASeq<-PrepareBam_EP(Samples=bam_file,
              SamplePath=PathToSamples,
            Ref_Transc="GTF",
            fileTransc=PathToGTF,
            cores=cores)


      print("Looking for AS events...")
      AllEvents_RNASeq<-EventDetection(SG_RNASeq, cores=cores, Path=TxtPath)

    },
    error=function(e){cat("\nERROR at sample",i, conditionMessage(e), "\n")},
    warning= function(w){cat("\nWARNING at sample ",i, conditionMessage(w), "\n")},
    finally= {print(cat("\n",samplelist[i], "\n"))}
    )
  }
}
