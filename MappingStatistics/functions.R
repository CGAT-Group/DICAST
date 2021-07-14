library(Rsamtools)
library(data.table)
library(GenomicFeatures)
library(tidyr)
library(GenomicAlignments)
library(pryr)

# extracts the exons by transcript of a .gtf file
readTranscripts <- function(path){
  
  print("Reading gtf file...")
  old <- Sys.time()
  
  txdb <- makeTxDbFromGFF(path)
  transcripts <- exonsBy(txdb,"tx",use.names=T)
  print(paste0(length(transcripts)," transcripts read."))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  return(transcripts)
}


# reads a fastq file and returns the headers as granges
fastqGRanges <- function(path){
  fastq <- readFastQHeadersAwk(path)
  
  print("Transforming to granges...")
  old <- Sys.time()
  
  granges <- GRanges(seqnames = Rle(fastq$transcript),
          ranges = IRanges(start=fastq$start,end=fastq$end),
          strand = Rle(strand(fastq$strand)),
          read_id = fastq$id,
          mate = fastq$mate,
          range_needs_check =fastq$range_needs_check
  )
  print(paste0(length(granges)," reads succesfully read..."))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  return(granges)
}

#reads a fastq files and returns formatted headers
readFastQHeaders <- function(path){
  
  #read in file
  print("Reading fastq file...")
  old <- Sys.time()
  
  fastq <-  scan(path, what="", sep="\n")
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  #extract header rows
  print("Get headers...")
  old <- Sys.time()
  
  header <-data.table(header = fastq[grep("@",fastq)])
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  #format
  print("Formatting...")
  old <- Sys.time()
  
  header[,strand:="*"]
  header$header <- gsub("@","",header$header)
  header <- separate(header,col=header,into=c("id","info"),sep= "\\/")
  header <- separate(header,col=info,into=c("transcript","mate1","mate2"),sep=";")
  header <- melt(data=header,id.vars=c("id","transcript","strand"),variable.name="mate",value.name="range")
  header$mate <- as.integer(gsub("mate","",header$mate))
  header$range <- sapply(header$range,function(x){
    strsplit(x,":")[[1]][2]
  })
  header <- separate(header,col=range,into=c("start","end"),sep="-")
  header[,start:=as.integer(start)]
  header[,end:=as.integer(end)]
  
  #deal with na ends (not correct currently)
  header$range_needs_check <-F
  header[is.na(header$end)]$range_needs_check<-T 
  header[is.na(header$end)]$end<-76
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  
  
  return(header)
}

#reads a .bam file
readBamFile <- function(bam_path){
  
  #read first mate
  print("reading first mate...")
  old <- Sys.time()
  param <- ScanBamParam( flag = scanBamFlag(isPaired = T,
                                            isProperPair = T,
                                            isUnmappedQuery = F,
                                            hasUnmappedMate = F,
                                            isSecondaryAlignment = F,
                                            isSupplementaryAlignment = F,
                                            isFirstMateRead = T
  ))
  mate_1 <- readGAlignments(bam_path, param=param, use.names = T)
  

  print("time:")
  print(Sys.time()-old)
  
  
  #read second mate
  print("reading second mate...")
  old <- Sys.time()
  param <- ScanBamParam( flag = scanBamFlag(isPaired = T,
                                            isProperPair = T,
                                            isUnmappedQuery = F,
                                            hasUnmappedMate = F,
                                            isSecondaryAlignment = F,
                                            isSupplementaryAlignment = F,
                                            isSecondMateRead = T
  ))
  mate_2 <- readGAlignments(bam_path,param=param, use.names = T)
  

  print("time:")
  print(Sys.time()-old)
  
  
  print("transforming to genomic ranges...")
  old <- Sys.time()
  #get the ranges
  bam <- list(mate_1=mate_1,mate_2=mate_2)

  bam_ranges_mate_1 <- granges(bam$mate_1)
  bam_ranges_mate_2 <- granges(bam$mate_2)
  

  bam_ranges_mate_1$mate <- 1
  bam_ranges_mate_2$mate <- 2

  

  print("time:")
  print(Sys.time()-old)
  
  
  print("parsing junctions from cigar...")
  old <- Sys.time()
  #get junctions from cigar

  bam_ranges_mate_1$junctions <-  extractAlignmentRangesOnReference(cigar(bam$mate_1),pos=start(bam$mate_1),drop.D.ranges = T)
  bam_ranges_mate_2$junctions <-  extractAlignmentRangesOnReference(cigar(bam$mate_2),pos=start(bam$mate_2),drop.D.ranges = T)
  #print(table(cigar(bam$mate_1)))

  print("time:")
  print(Sys.time()-old)
  
  
  # combine the mates into one granges object
  print("combining the two mates... ")
  old <- Sys.time()
  bam_ranges <- c(bam_ranges_mate_1,bam_ranges_mate_2)
  
  if(grepl("/", names(bam_ranges[1]), fixed = TRUE)){
    print("Parsing read names...")
    names(bam_ranges)<-unlist(strsplit(names(bam_ranges),"/",fixed=T))[ c(TRUE,FALSE) ]
  }
  if(grepl(" ",as.character(seqnames(bam_ranges[1])),fixed = TRUE)){
    print("Parsing chromosome names...")
    replace <- gsub("\\s.*","",as.character(seqnames(bam_ranges)))
    bam_ranges <- GRanges(seqnames = Rle(replace),
                          ranges = ranges(bam_ranges),
                          names= names(bam_ranges),
                          strand = strand(bam_ranges),
                          junctions = bam_ranges$junctions,
                          mate = bam_ranges$mate)
                    
  }
  bam_ranges$id <- names(bam_ranges)
  names(bam_ranges) <- paste(names(bam_ranges),bam_ranges$mate,sep="_")
  

  print("time:")
  print(Sys.time()-old)
  
  print("finished reading the bam file! Number of reads:")
  print(length(bam_ranges))
  return(bam_ranges)
}

readFastQHeadersAwk <- function(path){
  
  awk.script <- "read_fastq.awk"
  
  print("Reading fastq file with awk...")
  old <- Sys.time()
  
  header <- fread(text=system2("awk",args=c("-f",awk.script,path),stdout=TRUE))
  
  print("mem:")
  print(mem_used())
  print("time:")
  print(Sys.time()-old)
  return(header)
}

readConfig <- function(path){
  config <- fread(path,header = F)
  path_variables <- as.list(config$V2)
  names(path_variables)<- gsub(":","",config$V1)
  print("Path variables: ")
  print(path_variables)
  return(path_variables)
}

readArgs <- function(args){
  path_variables <- list(
    bam_path = args[1],
    fastq_path = args[2],
    gtf_path = args[3],
    genomic_reads_data = args[4],
    output_folder = args[5]
  )
  print("Path variables: ")
  print(path_variables)
  return(path_variables)
}

filterPairs <- function(bam_ranges){
  id_table <- data.table(id=bam_ranges$id,index=1:length(bam_ranges))
  counts <- id_table[,.(.N,index),by=id]
  double <- counts[N==2,index]
  old_length <- length(bam_ranges)
  new_length <- length(double)
  print(paste0(old_length-new_length," unpaired reads were removed"))
  return(bam_ranges[double])
}

removeImproperPairs <- function(bam_ranges){
  old_length <- length(bam_ranges)
  
  strand_table <- data.table(id=bam_ranges$id,index=1:length(bam_ranges),strand = as.character(strand(bam_ranges)))
  
  counts <- strand_table[,.N,by=c("id","strand")]
  same_strand <- unique(counts[N==2,id])
  strand_table <- strand_table[id %in% same_strand]
  bam_ranges <- bam_ranges[-strand_table$index]
  new_length <- length(bam_ranges)
  print(paste0(old_length-new_length," reads with mate on the same strand were removed"))
  
  chr_table <- data.table(id=bam_ranges$id,index=1:length(bam_ranges),chr = as.character(seqnames(bam_ranges)))
  counts <- chr_table[,.N,by=c("id","chr")]
  diff_chr <- unique(counts[N==1,id])
  chr_table <- chr_table[id %in% diff_chr]
  bam_ranges <- bam_ranges[-chr_table$index]
  print(paste0(new_length-length(bam_ranges)," reads with mate on different chromosomes were removed"))
  
  return(bam_ranges)
}

readBamMinimap <- function(bam_path,mate1_strand){
  old <- Sys.time()
  print("Reading bam file")
  param <- ScanBamParam(flag = scanBamFlag( isSecondaryAlignment = F,isSupplementaryAlignment = F))
  bam_ranges <- readGAlignments(bam_path, param = param, use.names = T)
  print(Sys.time()-old)
  
  print("parsing junctions from cigar...")
  old <- Sys.time()
  #get junctions from cigar
  junctions <-  extractAlignmentRangesOnReference(cigar(bam_ranges),pos=start(bam_ranges),drop.D.ranges = T)
  bam_ranges <- granges(bam_ranges)
  bam_ranges$junctions <- junctions
  junctions <- NULL
  print("time:")
  print(Sys.time()-old)
  
  old <- Sys.time()
  if(grepl("/", names(bam_ranges[1]), fixed = TRUE)){
    print("Parsing read names...")
    names(bam_ranges)<-unlist(strsplit(names(bam_ranges),"/",fixed=T))[ c(TRUE,FALSE) ]
  }
  if(grepl(" ",as.character(seqnames(bam_ranges[1])),fixed = TRUE)){
    print("Parsing chromosome names...")
    replace <- gsub("\\s.*","",as.character(seqnames(bam_ranges)))
    bam_ranges <- GRanges(seqnames = Rle(replace),
                          ranges = ranges(bam_ranges),
                          names= names(bam_ranges),
                          strand = strand(bam_ranges),
                          junctions = bam_ranges$junctions)
    
  }
  print("time:")
  print(Sys.time()-old)
  bam_ranges$id <- names(bam_ranges)
  
  
  old <- Sys.time()
  print("Filtering for pairs")
  bam_ranges <- filterPairs(bam_ranges)
  print("time:")
  print(Sys.time()-old)
  
  old <- Sys.time()
  print("Removing improper pairs")
  bam_ranges <- removeImproperPairs(bam_ranges)
  print("time:")
  print(Sys.time()-old)
  
  bam_ranges$mate <- 2
  bam_ranges[as.character(strand(bam_ranges))==mate1_strand]$mate <- 1
  names(bam_ranges) <- paste(names(bam_ranges),bam_ranges$mate,sep="_")
  
  print("finished reading the bam file! Number of reads:")
  print(length(bam_ranges))
  
  return(bam_ranges)
}

loadReferenceReads <-function(genomic_reads_data,fastq_path,gtf_path){
  if(!file.exists(genomic_reads_data)){
    #reading input files
    print("Generating genomic_reads...")
    reads <- fastqGRanges(fastq_path)
    transcripts <- readTranscripts(gtf_path)
    
    
    
    #correct false ranges
    print("Correcting ranges...")
    old <- Sys.time()
    
    read_width <- width(reads[Position(function(x) x==FALSE, reads$range_needs_check)])
    print(paste0("Read width: ",read_width))
    
    transcript_ids_for_range_parsing_1 <- as.character(seqnames(reads[reads$mate==1&reads$range_needs_check]))
    transcript_ids_for_range_parsing_2 <- as.character(seqnames(reads[reads$mate==2&reads$range_needs_check]))
    
    transcipts_for_range_parsing_1 <- transcripts[transcript_ids_for_range_parsing_1]
    transcipts_for_range_parsing_2 <- transcripts[transcript_ids_for_range_parsing_2]
    
    width_per_transcript_1 <- sum(width(transcipts_for_range_parsing_1))
    width_per_transcript_2 <- sum(width(transcipts_for_range_parsing_2))
    
    corrected_ranges_1 <- IRanges(start=1,end=pmin(width_per_transcript_1,read_width))
    corrected_ranges_2 <- IRanges(start=pmax(1,width_per_transcript_2-read_width+1),end=width_per_transcript_2)
    
    ranges(reads[reads$mate==1&reads$range_needs_check])<-corrected_ranges_1
    ranges(reads[reads$mate==2&reads$range_needs_check])<-corrected_ranges_2
    
    print(paste0(sum(reads$range_needs_check)," ranges corrected."))
    
    print(Sys.time()-old)
    
    
    
    
    #transform transcriptomic coordinates into genomic coordinates
    print("Parsing read coordinates to genomic...")
    old <- Sys.time()
    
    mapping <- mapFromTranscripts(reads,transcripts,ignore.strand=F)
    
    
    print(paste0(length(mapping)," read coordinates parsed..."))
    print("Unparsed reads:")
    print(reads[-mapping$xHits])
    print(Sys.time()-old)
    
    
    
    #get junctions
    print("Parsing junctions...")
    old <- Sys.time()
    
    mapping_list <- split(mapping,1:length(mapping))
    junctions <- intersect(mapping_list,transcripts[mapping$transcriptsHits],ignore.strand=T)
    
    print(Sys.time()-old)
    
    #generate genomic reads
    print("Generate genomic reads...")
    old <- Sys.time()
    
    genomic_reads <- GRanges(seqnames = seqnames(mapping),
                             ranges = ranges(mapping),
                             strand = strand(mapping),
                             mate = reads$mate
    )
    names(genomic_reads)<-paste(reads$read_id,reads$mate,sep="_")
    genomic_reads$junctions <- ranges(junctions)
    
    print(Sys.time()-old)
    
    #generate genomic reads
    print(paste0("Saving genomic reads as ", genomic_reads_data))
    old <- Sys.time()
    
    save(genomic_reads, file = genomic_reads_data)
    
    print(Sys.time()-old)
    return(genomic_reads)
  }else{
    
    print(paste0("Genomic reads are already computed. Loading data ",genomic_reads_data))
    old <- Sys.time()
    load(genomic_reads_data)
    print(paste0(length(genomic_reads)," reads loaded!"))
    print(Sys.time()-old)
    return(genomic_reads)
    
  }
}
