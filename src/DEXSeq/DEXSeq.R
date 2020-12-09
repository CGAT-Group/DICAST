suppressMessages(library(optparse))
suppressMessages(library(DEXSeq))
suppressMessages(library(BiocParallel))


#Parsing parameters
option_list = list(make_option(c("--gff"), type='character', default = NULL, help="the annotation gff file", metavar='character'),
		   make_option(c("--out"), type ='character',default = NULL, help="output directory",metavar='character'),
		   make_option(c("--ncores"), type ="integer", default=1,help="number of cores",metavar="integer"))



opt_parser = OptionParser(option_list=option_list) 
opt = parse_args(opt_parser)
BPPARAM = MulticoreParam(opt$ncores)


countFilesCase <- list.files(paste0(opt$out,"/tmp"),pattern="g1.txt$",full.names=T)
countFilesControl <- list.files(paste0(opt$out,"/tmp"),pattern="g2.txt$",full.names=T)


sampleData <- data.frame(row.names = c(basename(countFilesCase),basename(countFilesControl)),
			 condition = c(rep("case",length(countFilesCase)), rep("control",length(countFilesControl))))


#preparing dexseq object:
dxd <- DEXSeqDataSetFromHTSeq(c(countFilesCase,countFilesControl),sampleData=sampleData, design= ~ sample + exon + condition:exon,flattenedfile=opt$gff)

#testing for differential expression
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd,BPPARAM=BPPARAM)
dxd <- testForDEU(dxd,BPPARAM=BPPARAM)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition",BPPARAM=BPPARAM)

results <- DEXSeqResults(dxd)
write.table(results, paste0(opt$out, "/DEXSeqResults.txt"),sep="\t")

print("Finished DEXSeq analysis.")


