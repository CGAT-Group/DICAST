library(edgeR)
library(Rsubread)
library(optparse)

#The parameters parsing
option_list = list(make_option(c("--gtf"), type='character', default = NULL, help="the annotation gtf file", metavar='character'),
		   make_option(c("--casedir"), type='character', default = NULL, help="the dir of bam files", metavar='character'),
                   make_option(c("--controldir"), type='character', default = NULL, help="the dir of bam files", metavar='character'),
		   make_option(c("--output"), type='character', default = "output", help="the output suffix", metavar='character'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

setwd(opt$output)
if (is.null(opt$gtf)){
  print_help(opt_parser)
  stop("Please provide an annotation gtf file", call.=FALSE)
}

if (is.null(opt$case)){
  print_help(opt_parser)
  stop("Please provide case bam folder", call.=FALSE)
}

if (is.null(opt$control)){
  print_help(opt_parser)
  stop("Please provide control bam folder", call.=FALSE)
}

#Grouping bam files
cases <- list.files(opt$casedir,full.names=T, pattern="\\.bam$")
controls <- list.files(opt$controldir,full.names=T,pattern="\\.bam$")

bam_files = unlist(c(cases,controls))
Groups = as.integer(c(rep(1, length(cases)), rep (2,length(controls))))

#Reads counting 
fc <- featureCounts(bam_files, annot.ext=opt$gtf, isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_id", useMetaFeatures=FALSE, allowMultiOverlap=TRUE, isPairedEnd=TRUE)

#Preparing files for the analysis
y <- DGEList(counts=fc$count, group=Groups)
y$genes = fc$annotation

#Filtering
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

#Normalizing
y <- calcNormFactors(y)

#Estimating dispersion
design <- model.matrix(~Groups)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)

#Exon differential expression
fit <- glmQLFit(y, design, robust=TRUE)
qlf <- glmQLFTest(fit)

#Differential splicing analysis
sp <- diffSpliceDGE(fit, coef=2, geneid="GeneID", exonid="Start")

simes <- topSpliceDGE(sp, test="Simes", n=length(sp$exon.df.test))
ftest <- topSpliceDGE(sp, test="gene", n=length(sp$exon.df.test))
exon <- topSpliceDGE(sp, test="exon", n=length(sp$exon.df.test))
#Output
name1 = paste(opt$output, "edgeR_simes.txt", sep="/")
name2 = paste(opt$output, "edgeR_ftest.txt", sep="/")
name3 = paste(opt$output, "edgeR_exon.txt", sep="/")

write.table(simes, sep='\t', file=name1)
write.table(ftest, sep='\t', file=name2)
write.table(exon, sep='\t', file=name3)
