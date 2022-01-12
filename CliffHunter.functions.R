################################################
#      Functions to MPNSTs samples analysis    #
################################################


source(file = "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/utils.R" )

#getSVposition

getSVposition <- function(bampath, ref_genome = "hg19", chr.style = "Ensembl", 
                          tilewidth = 100e6, percent.value=0.2, diff.value=30, min.mapqual = 30, verbose = TRUE){
  
  now.msg("Getting breakpoint positions...", verbose = verbose)
  
  msg("Creating genomic windows", verbose = verbose)
  genome<- filterChromosomes(getGenome(ref_genome))
  seqlevelsStyle(genome) <- chr.style
  
  #length of chr
  chr.lens<- setNames(end(genome), seqnames(genome))
  
  #create the windows
  windows <- tileGenome(seqlengths = chr.lens, tilewidth= tilewidth, cut.last.tile.in.chrom = TRUE)
  
  msg("Computing coverage differences. This may take a while...", verbose = verbose)
  
  sv.position <- data.frame()
  for(i in 1:length(windows)){
    msg(paste0("Exploring window ", i, " out of ", length(windows)), verbose = verbose)
    
    #get BamCoverage as a vector selecting with window_num
    coverage <- bamCoverage(bampath = bampath,windows[i], mapqual = min.mapqual ,verbose = FALSE)[1]
    
    #Compute the diference in coverage between each position and the one before
    difs <- (coverage-c(coverage[1], coverage[1:length(coverage)-1]))
    #detect larges differences, cliffs
    cliffs <- which(abs(difs)>diff.value)
    
    difs <- difs[cliffs]
    
    #compute the percentage of diffrence relative to the cliff top
    #IMPORTANT: for positive differences, the base of the cliff is one base before the detected difference, and so, the cliff  top is at the difference position
    #For negative differences, the cliff base is at the detected difference and the top is one base before
    percentage <- numeric(length(cliffs))
    
    percentage[difs>0] <- difs[difs>0]/coverage[cliffs[difs>0]]
    cov.base.before <- c(coverage[1], coverage[1:length(coverage)-1])
    percentage[difs<0] <- abs(difs[difs<0]/cov.base.before[cliffs[difs<0]])
    
    #Filter cliffs by minimum percentage
    valid.percentage <- which(percentage>percent.value)
    cliffs <- cliffs[valid.percentage]
    difs <- difs[valid.percentage]
    percentage <- percentage[valid.percentage]
    
    #Since the positive cliffs are detected on the top of the cliff, remove one position to get the base of the cliff
    cliffs[difs<0] <- cliffs[difs<0]-1
    
    #cliffs direction, if the cliff is forward(from the left)  or reverse(from the rigth) 
    direction <- character(length(difs))
    direction[difs>0]<-"rv"
    direction[difs<0]<-"fwd"
    
    #Transform to genomic position
    cliffs <- cliffs + start(windows[i]) -1
    
    sv.position <- rbind(sv.position, 
                         data.frame(chr=rep(as.character(seqnames(windows[i])),length(cliffs)),
                                    start=cliffs,
                                    end=cliffs,
                                    difference=difs,
                                    percentage=percentage,
                                    direction=direction)
    )
    
    
  } 
  
  
  now.msg("Breakpoint positions obtained", verbose = verbose)
  return(toGRanges(sv.position))
}


#annotation cliffs brakepoints

bpAnnotation <- function(sv_data, org.db= org.Hs.eg.db, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,txdb_columns="SYMBOL", txdb_keytype="ENTREZID"){
  probes_int <-  filterChromosomes(sv_data)
 
  loc <- locateVariants(query=probes_int, subject=txdb, AllVariants())
  names(loc) <- seq(1:length(loc))
  annot <- toDataframe(loc)
  
  for(i in 1:length(annot$start)){
    n <- annot$start[i]
    if(n %in% start(sv_data)){
      annot$difference[i] <- sv_data$difference[which(n==start(sv_data))]
      annot$percentage[i] <- sv_data$percentage[which(n==start(sv_data))]
      annot$direction[i] <- as.character(sv_data$direction[which(n==start(sv_data))])
    } 
  }

  annot<-annot[,c(1,2,3,13,14,4,10,15)]
  data_annot <- unique(cbind(annot,select(org.db, annot$GENEID, columns = txdb_columns, keytype = txdb_keytype)))
  data_annot <- unique(data_annot[,c("chr","start","end","SYMBOL","percentage","difference","direction","LOCATION")])
  data_annot <- toGRanges(data_annot)
  return(data_annot)
  
}

#########################
#   Soft clipped reads  #
#########################

#we need to filter cliffs where soft clipped sequences are

getCliffSoftSeq<- function(cliffs,
                           bam.file,
                           percent.nucleotide=0.8,
                           bam.fields = c("pos","seq"),
                           min.cons.numseqs = 4,
                           verbose = TRUE){
  
  now.msg("Starting getCliffSoftSeq function...", verbose = verbose)
  
  msg("Filtering bamfile...", verbose = verbose)
  cliffs$clipped.seq <- character(length(cliffs))
  
  #parameter selection of the bamfile from cliffs_gr output
  param <- ScanBamParam(which = cliffs, what = bam.fields)
  bamdata <- readGAlignments(bam.file, param = param)
  
  #we delete those reads with no soft clipping or soft clipping in both sides of it
  bamdata <- bamdata[which((nchar(cigar(bamdata)) - nchar(gsub("S", "", cigar(bamdata)))==1)),]
  cliffs <- cliffs[start(cliffs)%in% start(bamdata) | start(cliffs) %in% end(bamdata),]
  
  
  msg("Bamfile filtered...", verbose = verbose)
  
  # library(profvis)
  # profvis({
  
  now.msg("Extracting soft clipped sequences...", verbose = verbose)
  for(i in seq_len(length(cliffs))){

    s <-start(cliffs)[i]
    bam.soft.seq <- bamdata[which(s == start(bamdata)|s == end(bamdata)),]
    
    #we select only those reads which have been clipped in the desired cliff.    
    cliff.orientation <- as.character(cliffs$direction[i])
    #as.character(cliffs$cliff.direc[which(start(cliffs) == s)])
    
    #vec[i] <- c(which(sapply(cigar(zy), function(x) hasLetterAt(x,"S",nchar(x),fixed=T)) & cliff.orientation =="fwd"),which(!(sapply(cigar(zy), function(x) hasLetterAt(x,"S",nchar(x),fixed=T))) & cliff.orientation=="rv"))
    
    if(cliff.orientation =="fwd"){
      bam.soft.seq <- bam.soft.seq[grepl("S$", cigar(bam.soft.seq)),]
      
    }else{
      bam.soft.seq<-bam.soft.seq[!(grepl("S$", cigar(bam.soft.seq))),]
      #zy<-zy[!(sapply(cigar(zy), function(x) hasLetterAt(x,"S",nchar(x),fixed=T))),]
    }
    
    if(length(bam.soft.seq)==0) next
    
    #vector with all softclipping numbers and selection of the highest value
    n <-as.numeric(gsub('.*?([[:digit:]]+)S ?.*','\\1', cigar(bam.soft.seq)))
    
    
    #we can have different reads with similar clipped sequences of the same maximal length.That's why we decided to select the first one
    ###ATENTION### We will modify this method and we will select the longest consensus sequence.
    if(length(n)>=5){
      if(cliff.orientation =="fwd"){
        # soft.clip.seq<-substr(mcols(bam.soft.seq)$seq[maxpos], nchar(mcols(bam.soft.seq)$seq[maxpos])-(max(n)[1])+1, nchar(mcols(bam.soft.seq)$seq[maxpos]))
        soft.clip.seq<-substr(mcols(bam.soft.seq)$seq, nchar(mcols(bam.soft.seq)$seq)-(n)+1, nchar(mcols(bam.soft.seq)$seq))
      }else{
        #soft.clip.seq<-substr(mcols(bam.soft.seq)$seq[maxpos], 1, max(n)[1])
        soft.clip.seq<-substr(mcols(bam.soft.seq)$seq, 1, n) 
      }
      
      #consensus seq
      dna <- DNAStringSet(soft.clip.seq)
      DNA <- AlignSeqs(dna,refinements = 1, verbose=F)
      seq.align <-t(data.frame(strsplit(as.character(DNA),"")))
      #consensus<- as.character()
      
      msg(paste("Getting consensus sequence", i,"of",length(cliffs),"..."), verbose = verbose)
      cons.nuc <- character()
      cons.pct <- numeric()
      cons.numseqs <- numeric()
      
      for(j in seq_len(ncol(seq.align))){
        nuc <- seq.align[,j]
        
        nucleotide<- sort(table(nuc[!(nuc=="-")]), decreasing = TRUE)#[1]
        cons.numseqs <- c(cons.numseqs, sum(nucleotide))
        
        #we calculate the percentage of the presence of each nucleotide
        #in one specific position
        nuc.percent<-nucleotide/sum(nucleotide)
        cons.pct <- c(cons.pct, nuc.percent[1])
        
        cons.nuc <- c(cons.nuc, nucleotide[1])
        
      }
      
      #we filter consensus seq
      if(cliff.orientation =="fwd"){
        seq.start <- which(cons.numseqs>min.cons.numseqs)[1]
        seq.start <-na.omit(seq.start)
        seq.end <- length(cons.numseqs)
      }else{
        seq.start <- 1
        seq.end <- which(cons.numseqs>min.cons.numseqs)[length(which(cons.numseqs>min.cons.numseqs))]
      }
      
      if(length(seq.end) ==0 |length(seq.start)==0)next
      
      if(seq.end>seq.start) {
        if(seq.end-seq.start<10) next 
        if(length(which(cons.pct[seq.start:seq.end]>0.9))<percent.nucleotide*(seq.end-seq.start)) next 
        cliffs$clipped.seq[i] <- paste0(names(cons.nuc[seq.start:seq.end]), collapse="")
      } else {
        cliffs$clipped.seq[i]<-NA
      }
    } else{
      cliffs$clipped.seq[i]<-NA
    }
    
  }
  
  cliffs<- cliffs[which(cliffs$clipped.seq!=""),]
  msg("Soft clipped sequences extracted", verbose = verbose)
  now.msg("getCliffSoftSeq done", verbose = verbose)
  
  return(cliffs)
}


##########################################
# Soft clipped sequences alignment: BWA  #
##########################################

#Once we have our soft clipped sequences, we proceed to know whether they matchor not in different parts of the reference genome,
#therefore they are insertions|translocations, or whether they are deletions

#First we save a Fastq file with Cliff positions as name of the soft sequence obtained
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

#We align with BWA
fastqToBAM <- function(fastq.file, threads.num=2, ref.genome) {
  
  #Map Them against the genome
  
  bwa.cmd <- paste0("bwa mem -t ", threads.num, " ", ref.genome, " ", fastq.file, " > ", paste0(fastq.file, ".sam"))
  message("Mapping... (", bwa.cmd, ")")
  system(command=bwa.cmd, wait=TRUE)
  
  #And convert the BAM to SAM
  message("Converting...")
  bam.cmd <- paste0("samtools view -Sb ", paste0(fastq.file, ".sam"), " > ", paste0(fastq.file, ".bam"))
  system(command=bam.cmd, wait=TRUE)
  message("Sorting...")
  sort.cmd <- paste0("samtools sort ", paste0(fastq.file, ".bam"), " > ", paste0(fastq.file, ".bam.sorted.bam"))
  system(command=sort.cmd, wait=TRUE)
  file.remove(paste0(fastq.file, ".bam"))
  file.rename(from=paste0(fastq.file, ".bam.sorted.bam"), to=paste0(fastq.file, ".bam"))
  index.cmd <- paste0("samtools index ",  paste0(fastq.file, ".bam"))
  message("Indexing...")
  system(command=index.cmd, wait=TRUE)  
}