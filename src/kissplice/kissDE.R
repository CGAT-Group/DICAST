#!/usr/bin/Rscript

suppressMessages(library(kissDE))
suppressMessages(library(optparse))

#Parse arguments
option_list = list(make_option(c("--counts"), type='character', default = NULL, help="the output file of kissplice2refgenome (full path)", metavar='character'),
                   make_option(c("--cores"), type='integer', default = 1, help="the number of cores"),
                   make_option(c("--caseprefix"),type='character', default=NULL,help="prefix of all files contained in casefolder",metavar='character'),
		   make_option(c("--controlprefix"),type='character', default=NULL,help="prefix of all files contained on controlfolder",metavar='character'),
                   make_option(c("--casefastq"),type='character', default=NULL,help="full path to casefastq folder",metavar='character'),
                   make_option(c("--controlfastq"),type='character', default=NULL,help="full path to controlfastq folder",metavar='character'),
                   make_option(c("--out"), type='character',default=NULL,help="the full path to the output folder",metavar='character'))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#handle empty count file
if(is.null(opt$counts)) stop(paste0("the handed count file ",opt$counts," does not exist or is empty.\n Stopping DS analysis."))
cores<- opt$cores

#load count data for AS events from kissplice2refgenome output
countsData<-kissplice2counts(opt$counts,k2rg=T)

#get number of case/control files
ncasefiles<-length(list.files(opt$casefastq))
ncontrolfiles<-length(list.files(opt$controlfastq))
conditions <- c(rep(opt$caseprefix,ncasefiles),rep(opt$controlprefix,ncontrolfiles))


results <- diffExpressedVariants(countsData, conditions,nbCore=cores)

writeOutputKissDE(results,adjPvalMax=.05, dPSImin=.1,output=paste0(opt$out,"/kissDE-output.tab"))



