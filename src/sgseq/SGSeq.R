suppressMessages(library(optparse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(SGSeq))

#Parsing parameters
option_list = list(make_option(c("--gtf"), type='character', default = NULL, help="the annotation gtf file", metavar='character'),
                   make_option(c("--cores"), type='integer', default = 1, help="the number of cores"),
                   make_option(c("--path_to_bam"), type='character', default = NULL, help="the path to bam file", metavar='character'),
                   make_option(c("--out"), type='character',default=NULL,help="the full path to the output folder",metavar='character'))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide an annotation gtf file", call.=FALSE)
}

if (is.null(opt$path_to_bam)){
  print_help(opt_parser)
  stop("Please provide a path to bam files", call.=FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Please provide a full path to the output folder", call.=FALSE)
}
cores=opt$cores


#Building the config dataframe
sample_name <- basename(opt$path_to_bam)
config <- data.frame(sample_name=sample_name, file_bam=opt$path_to_bam)
si <- getBamInfo(sample_info=config, yieldSize=NULL,cores=cores)

#Reading the custom annotation
txdb <- makeTxDbFromGFF(opt$gtf, format="gtf")
txf_ucsc <- convertToTxFeatures(txdb)
sgf_ucsc <- convertToSGFeatures(txf_ucsc)

#Annotating splice variants from the custom annotation
sgfc_ucsc <- analyzeFeatures(si, features = txf_ucsc, cores=cores)
sgvc_ucsc <- analyzeVariants(sgfc_ucsc, cores=cores)
sgvc_ucsc_fpkm <- variantFreq(sgvc_ucsc) 
result <- cbind(data.frame(mcols(sgvc_ucsc)), data.frame(sgvc_ucsc_fpkm))
write.csv(result, paste0(opt$out,"/SGSeq_from_gtf.csv"))

#Annotating splice variants de novo
sgfc_pred <- analyzeFeatures(si, cores=cores)
sgfc_pred <- annotate(sgfc_pred, txf_ucsc)
sgvc_pred <- analyzeVariants(sgfc_pred, cores=cores)
sgvc_pred_fpkm <- variantFreq(sgvc_pred) 
result <- cbind(data.frame(mcols(sgvc_pred)), data.frame(sgvc_pred_fpkm))
write.csv(result, paste0(opt$out,"/SGSeq_denovo.csv"))

print("----------------- SGSeq-script finished --------------------")
