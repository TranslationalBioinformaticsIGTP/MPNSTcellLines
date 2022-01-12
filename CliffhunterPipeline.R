###################################################
#       SV analysis using CliffHunteR             #
###################################################
# Packages Needed
library(CopyNumberPlots)
library(bamsignals)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(GenomicAlignments)
library(DECIPHER)
library(AnnotationHub)
library(StructuralVariantAnnotation)
############## Functions ################
source(file = "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/CliffHunter.functions.R")

#.unlist()
.unlist <- function (x){
  ## do.call(c, ...) coerces factor to integer, which is undesired
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#Gene annotation 
getGeneAnnot <- function(df, txdb ){
  loc <- VariantAnnotation::locateVariants(query = toGRanges(df),subject = txdb,AllVariants())
  k <- loc$GENEID
  tx2gene <- AnnotationDbi::select(x = orgdb,
                                   keys = k,
                                   columns = "SYMBOL",
                                   keytype = "ENTREZID")
  tx2genes <- unique(tx2gene)
  loc$SYMBOL <-""
  
  for(l in seq_len(length(unique(loc$GENEID)))){
    gid <- unique(loc$GENEID)[l]
    gid.pos <- which(gid== loc$GENEID)
    if(is.na(gid)){
      loc$SYMBOL[gid.pos] <- NA
    } else{
      loc$SYMBOL[gid.pos] <- tx2genes$SYMBOL[which(gid ==tx2genes$ENTREZID)]
      
    }
    
  }
  return(loc)
}

#################### Parameters & Directories #####################
# cell.lines
sample.names <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/Sample.info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample.names <- sample.names[,2]
cell.line <- "S462"
# cell.line <- "NMS.2"
# cell.line <- "ST88-14"

ref.genome.fq <- "/imppc/labs/eslab/mmagallon/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref.genome <- "hg38"
orgdb <- org.Hs.eg.db

#repeatMasker 
rptmskr <- read.table(file.path("/imppc/labs/eslab/mmagallon/Annotatiions/RepeatMasker/rmskr_hg38.bed"),header = FALSE, stringsAsFactors = FALSE)
rptmskr <- toGRanges(rptmskr, genome = ref.genome)
rptmskr <-filterChromosomes(rptmskr)

#CNV UCSC track
cnv.ucsc.hg38 <- read.table("/imppc/labs/eslab/mmagallon/Annotatiions/DGV_population.bed", sep = "\t",header = FALSE)
cnv.ucsc.hg38 <- toGRanges(cnv.ucsc.hg38)
unique(cnv.ucsc.hg38$V9)
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "139,69,19"] <- "Loss_Gain"
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "200,0,0"] <- "Loss"
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "200,0,200"] <- "Inv"
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "0,0,200"] <- "Gain"
cnv.ucsc.hg38 <- sortSeqlevels(cnv.ucsc.hg38)
cnv.ucsc.hg38 <- sort(cnv.ucsc.hg38)

#Selection of interesting genes
gene.markers <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_cellLines/WGS","special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
gene.markers.gf <- toGRanges(gene.markers$V2)

mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
#seleccted Genes
gene.markers.gf <- gene.markers.gf[gene.markers.gf$Genes %in% c("NF1","SUZ12","EED",
                                                                "EZH1",
                                                                "EZH2",
                                                                "CDKN2A",
                                                                "TP53",
                                                                "PTEN",
                                                                "RB1",
                                                                "ERBB2",
                                                                "ERBB3",
                                                                "KIT",
                                                                "HGF",
                                                                "EGFR",
                                                                "MET",
                                                                "SOX9",
                                                                "BIRC5",
                                                                "PDGFRA",
                                                                "AURKA",
                                                                "AURKB",
                                                                "NF2",
                                                                "LZTR1",
                                                                "CXCR4",
                                                                "TWIST1",
                                                                "BCR")]
regions.selected <- gene.markers.gf + 1e6
regions.df <- toDataframe(regions.selected)

# Directories
execution.dir <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS"

#### Execution CliffHunter ####
sample.names <- sample.names[-1]
for(i in seq_len(length(sample.names))){
  cell.line <- sample.names[i]
  message("CliffHunteR of ", cell.line, " cell-line")
  #Directoreies
  result.dir <- file.path(execution.dir,"results",cell.line)
  cliff.dir <- file.path(result.dir,"CliffHunteR")
  if(!file.exists(cliff.dir)) dir.create(cliff.dir)
  complete.bampath <- file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/results", cell.line,"BAM", paste0(cell.line, ".bam"))
  
  
  bampath <- file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/results/", cell.line,"CliffHunteR", paste0(cell.line, "_selected_regions.bam"))
  
  if(!file.exists(bampath)){
    
    write.table(regions.df[,c(1:3)], file.path(cliff.dir, "selected.regions.bed"),quote = FALSE,col.names = FALSE, row.names = FALSE,sep = "\t")
    
    # Selecting the specifcic regions to look for SV
    cmd <- paste0("samtools view ", " -b -L ", cliff.dir, "/selected.regions.bed ", complete.bampath, " > ", cliff.dir,"/", cell.line, "_selected_regions.bam")
    system(cmd)
    
    #Indexing the bam
    cmd <- paste0("samtools index ", cliff.dir,"/", cell.line, "_selected_regions.bam")
    system(cmd)
    
  }
  
  ################### Cliff analysis ###############
  bampath <- file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/results/", cell.line,"CliffHunteR", paste0(cell.line, "_selected_regions.bam"))
  
  SV <- getSVposition(bampath = bampath, ref_genome = ref.genome,
                      chr.style = "UCSC",
                      min.mapqual =10,
                      diff.value =4,
                      percent.value = 0.25)
  names(SV) <- paste0(as.character(seqnames(SV)),":", start(SV),"-",end(SV))
  
  #Annotation Cliffs
  SV.annot <- bpAnnotation(sv_data = SV, org.db =  org.Hs.eg.db, txdb = TxDb.Hsapiens.UCSC.hg38.knownGene)
  names(SV.annot) <- paste0(as.character(seqnames(SV.annot)),":", start(SV.annot),"-",end(SV.annot))
  
  #Filttering common CNVs annotated in UCSC
  sv<- SV.annot
  #find overlaps to annotate
  mcols(sv)$DGV_symbol <-""
  mcols(sv)$DGV_SV <-""
  str.end <- c("start", "end")
  for(i in seq_len(length(str.end))){
    fnd <- findOverlaps(query = cnv.ucsc.hg38,sv,type = str.end[i])
    fnd <- data.frame(fnd)
    
    sv[fnd$subjectHits]$DGV_symbol <- cnv.ucsc.hg38[fnd$queryHits]$V4
    sv[fnd$subjectHits]$DGV_SV <- cnv.ucsc.hg38[fnd$queryHits]$V9
  }
  
  
  #RepeatMaskerAnnotation
  sv$repeatMasker <- ""
  
  fnd <- findOverlaps(query = rptmskr,sv,type = "any")
  fnd <- data.frame(fnd)
  sv[fnd$subjectHits]$repeatMasker <- rptmskr[fnd$queryHits]$V4
  data.frame(sv)
  
  # #Filter those SV commonly found in cnv.ucsc.hg38
  # sv.filt <- sv[sv$DGV_symbol ==""]
  # sv.filt <- sv.filt[sv.filt$repeatMasker != "(AT)n"]
  SV <- getCliffSoftSeq(cliffs = sv, bam.file = bampath,
                        min.cons.numseqs =4,
                        percent.nucleotide = 0.8 )
  
  
  
  
  # Getting Cliffs Sequences: This funtion filter and select the best Sequences.
  # That is why we loss some of the interesting cliffs.
  # SV <- getCliffSoftSeq(cliffs = SV, bam.file = bampath,
  #                       min.cons.numseqs = 2,
  #                       percent.nucleotide = 0.8 )
  
  # Writting a fastq file containing all softclipped reads detected
  sv.data <- data.frame(name = names(SV), seq = SV$clipped.seq)
  writeFasta(data = sv.data, filename = file.path(cliff.dir,paste0("soft_seq_0.2_0.8_", cell.line,".fq")))
  
  #Aligment of softclipped sequences selected
  fastqToBAM( fastq.file = file.path(cliff.dir, paste0("soft_seq_0.2_0.8_", cell.line,".fq")), ref.genome = ref.genome.fq, threads.num = 2 )
  
  #Read bam file generated to associate the breakends.
  new.bp.bam <- scanBam(file = file.path(cliff.dir, paste0("soft_seq_0.2_0.8_", cell.line,".fq.bam")))
  bam.field <- names(new.bp.bam[[1]])
  
  list <- lapply(bam.field, function(y) .unlist(lapply(new.bp.bam, "[[", y)))
  bam.df <- do.call("data.frame", list)
  head(bam.df)
  names(bam.df) <- bam.field
  # chr.9.test <- bam.df[grepl(pattern = "chr9",x = bam.df$qname),]
  # chr.9.test.gr <- toGRanges(chr.9.test$qname)
  # mcols(chr.9.test.gr) <- chr.9.test[,2:ncol(chr.9.test)]
  # chr.9.test.gr <- sort(chr.9.test.gr)
  
  all.gr <- toGRanges(bam.df$qname)
  mcols(all.gr) <- bam.df[,2:ncol(bam.df)]
  all.gr <- sortSeqlevels(all.gr)
  all.gr <- sort(all.gr)
  names(all.gr) <- paste0(as.character(seqnames(all.gr)),":", start(all.gr), "-", end(all.gr))
  
  
  mate.pos <- paste0(all.gr$rname,":",all.gr$pos, "-", all.gr$pos)
  mate.pos[grepl("NA",mate.pos)]<- ""
  names(mate.pos) <- names(all.gr)
  
  SV$mate.pos <- ""
  sv.names <- names(SV)[names(SV)%in% names(mate.pos)]
  mate.names <- names(mate.pos)[names(mate.pos)%in%names(SV)]
  names(mate.names)<- mate.names
  mate.names <- mate.names[sv.names]
  mate.pos <- mate.pos[mate.names]
  SV[sv.names]$mate.pos<- mate.pos
  SV <- unique(SV)
  
  #Rename of genes based on our special.genes
  gene.markers <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_cellLines/WGS/","special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
  gene.markers.gf <- toGRanges(gene.markers$V2)
  mcols(gene.markers.gf) <- gene.markers$V1
  colnames(mcols(gene.markers.gf)) <- "Genes"
  gene.markers.gf <- sort(gene.markers.gf)
  
  for(l in seq_len(length(SV))){
    ss <- SV[l]
    s <- subsetByOverlaps(ss,gene.markers.gf)
    if(length(s)!=0){
      g <- subsetByOverlaps(gene.markers.gf,s)
      g <- g$Genes
      SV$SYMBOL[l] <- g
    }
  }
  mcols(SV) <- mcols(SV)[c(1:8,10,9)]
  
  #saving the results
  # save(SV,file = file.path(cliff.dir,paste0(cell.line,"_CliffHunter_CNV_0.6_4.RData")))
  save(SV,file = file.path(cliff.dir,paste0(cell.line,"_CliffHunter_CNV_0.25_4.RData")))
  
  SV.df <- toDataframe(SV)
  # write.table(SV.df,file = file.path(cliff.dir,paste0(cell.line,"_CliffHunter_CNV_0.6_4.csv")),sep = "\t", col.names = TRUE, row.names = TRUE)
  write.table(SV.df,file = file.path(cliff.dir,paste0(cell.line,"_CliffHunter_CNV_0.25_4.csv")),sep = "\t", col.names = TRUE, row.names = TRUE)
  
}
##### RData Filtering Specific Genes  #####
i = 1
SV.selregion.list <- list()
for(i in seq_len(length(sample.names))){
  cell.line <- sample.names[i]
  
  #Directoreies
  result.dir <- file.path(execution.dir,"results",cell.line)
  cliff.dir <- file.path(result.dir,"CliffHunteR")
  
  # load(file.path(cliff.dir, paste0(cell.line,"_CliffHunter_CNV_0.6_4.RData")))
  load(file.path(cliff.dir, paste0(cell.line,"_CliffHunter_CNV_0.25_4.RData")))
  
  #Selecting only interesting gene regions
  
  SV.selregion <- subsetByOverlaps(SV,gene.markers.gf)
  SV.selregion.list[[cell.line]] <- SV.selregion
  SV.df <- toDataframe(SV.selregion)
  write.table(SV.df,file = file.path(cliff.dir,paste0(cell.line,"_CliffHunter_CNV_0.25_4_Selected_Regions.csv")),sep = "\t", col.names = TRUE, row.names = TRUE)
  
}

ss <- unlist(as(SV.selregion.list, "GRangesList"))
names(ss) <- paste0(as.character(seqnames(ss)), ":", start(ss),"-", end(ss))
dup.svs <- names(ss)[which(duplicated(ss))]
i=1

for(i in seq_len(length(SV.selregion.list))){
  cell.line <- names(SV.selregion.list)[i]
  
  #Directoreies
  result.dir <- file.path(execution.dir,"results",cell.line)
  cliff.dir <- file.path(result.dir,"CliffHunteR")
  
  sv.data <- SV.selregion.list[[cell.line]]
  length(sv.data)
  sv.data <- sv.data[!names(sv.data) %in% dup.svs]
  length(sv.data)
  sv.data <- toDataframe(sv.data)
  
  #Annotation of SV
  loc <- getGeneAnnot(df = sv.data, txdb = txdb)
  sv.data$Genes <- ""
  for(l in seq_len(length(unique(names(loc))))){
    n.loc<- unique(names(loc))[l]
    sv.data[rownames(sv.data) == n.loc,"Genes"] <- paste0(unique(loc$SYMBOL[names(loc) == n.loc]),collapse =",")
    
  }
  # Annot mate possition
  if(length(sv.data$mate.pos[sv.data$mate.pos !=""])>0){
    
  }
  
  if(length(sv.data$mate.pos[sv.data$mate.pos !=""])>0){
    mate.pos.annot <-toGRanges(sv.data$mate.pos[sv.data$mate.pos !=""])
    
    names(mate.pos.annot) <- sv.data$mate.pos[sv.data$mate.pos !=""]
    
    loc <- getGeneAnnot(df = toDataframe(mate.pos.annot), txdb = txdb)
    sv.data$mate.pos.Genes <- ""
    
    for(l in seq_len(length(unique(names(loc))))){
      n.loc<- unique(names(loc))[l]
      sv.data[sv.data$mate.pos%in%n.loc,"mate.pos.Genes"] <- paste0(unique(loc$SYMBOL[names(loc) == n.loc]),collapse =",")
      
    }
  }
  
  
  
  write.table(sv.data,file = file.path(cliff.dir,paste0(cell.line,"_CliffHunter_CNV_0.25_4_Selected_Regions_FilteredCommonSVinSamples.csv")),sep = "\t", col.names = TRUE, row.names = TRUE)
  
}


