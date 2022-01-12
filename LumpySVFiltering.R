########################################
#           Lumpy Filtering            #
########################################
# Packages Needed
library(CopyNumberPlots)
library(CNVfilteR)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
# options(timeout = 300)
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
###################### Functions ############################
source(file = "/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/loadingLumpySVs.R")

###################### Parameters ############################
# Directories
execution.dir <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS"
# gene.markers <- toGRanges(data.frame(chr= c("chr9", "chr17", "chr17", "chr17", "chr11"), start = c(21967753,31094927, 31937007,7668230, 86245050), end = c(21975098, 31377677, 32001038, 7687366,86278810), symbol = c("CDKN2A", "NF1", "SUZ12", "TP53", "EED")))

#Cell lines
sample.names <- read.table(file.path(execution.dir,"Sample.info.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sample.names <- sample.names$Sample.name

#Genome
ref.genome.fq <- "/imppc/labs/eslab/mmagallon/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
genome <- "hg38"
orgdb <- org.Hs.eg.db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

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
gene.markers.gf <- toGRanges(gene.markers$V2,genome = genome)
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6
regions.df <- toDataframe(regions.selected)
gene.markers.gf <- gene.markers.gf[gene.markers.gf$Genes !="BCR"]

# TSG  COSMIC list
all.genes.cosmic <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS","Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]

all.genes.cosmic[545,]
tumor.sup.gr <- all.genes.cosmic.gr[grepl("TSG", all.genes.cosmic.gr$Role.in.Cancer)]
# tumor.sup.gr <- toGRanges(tumor.sup$Genome.Location)
# mcols(tumor.sup.gr) <- tumor.sup[,c(1,2,5,6:20)]
length(tumor.sup.gr)
# gene.markers.gf <- tumor.sup.gr

# gene.markers.gf <-all.genes.cosmic.gr
#Fusion Genes Cosmic
read.csv
table(all.fusion.genes.cosmic$PRIMARY_HISTOLOGY)
all.fusion.genes.cosmic <- read.csv(file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS","CosmicFusionExport.tsv"), sep = "\t",header = T,comment.char = "")

all.fusion.genes.cosmic <- data.frame(chr = paste0("chr",all.fusion.genes.cosmic$X3._CHROMOSOME),
                                      start = all.fusion.genes.cosmic$X3._GENOME_START_FROM,
                                      end = all.fusion.genes.cosmic$X3._GENOME_STOP_FROM,
                                      gene = all.fusion.genes.cosmic$X3._GENE_NAME,
                                      fusion.gene = all.fusion.genes.cosmic$X5._GENE_NAME,
                                      fusion = all.fusion.genes.cosmic$TRANSLOCATION_NAME,
                                      histology = all.fusion.genes.cosmic$PRIMARY_HISTOLOGY)
all.fusion.genes.cosmic <- all.fusion.genes.cosmic[all.fusion.genes.cosmic$chr !="chrNA",]
all.fusion.genes.cosmic.gr <- toGRanges(all.fusion.genes.cosmic, genome = "hg38")



######################## Directories #####################################
for(sample in sample.names){
  # images.dir <- file.path(execution.dir,"results", sample, "CNVkit/ExcludedRegions2","images")
  lumpy.dir <- file.path(execution.dir,"results", sample, "Lumpy")
  if(!file.exists(lumpy.dir))dir.create(lumpy.dir)
  # if(!file.exists(images.dir))dir.create(images.dir)
}

###########  Lumpy file ################
vcf.file <-file.path(execution.dir,"results/Lumpy_smoove_hg38/MPNST_cellLines_mw6-smoove.genotyped.vcf.gz")
genome <- "hg38"

all.sv.samples <- list()
for(i in seq_len(length(sample.names))){
  sample.name <- sample.names[i]
  message("SV from ", sample.name)
  svs <-loadingLumpySVs(vcf.file = vcf.file, genome= genome, sv.type = "all", sample.name = sample.name)
  
  svs$ICR$SVTYPE <- "ICR"
  svs$INV$SVTYPE <- "INV"
  svs <-  unlist(as(svs, "GRangesList"))
  names(svs) <- paste0(seqnames(svs),":",start(svs),"-",end(svs))
  # svs.nm <- strsplit(names(svs),"\\.")
  # svs.nm <- unlist(lapply(svs.nm, function(x){ x[2]}))
  # names(svs) <-  svs.nm
  all.sv.samples[[sample.name]] <- sort(svs)
}

# Getting the translocation of the region
transloc <- GRanges()
inv <- GRanges()
dup <- GRanges()
del <- GRanges()

for(i in seq_len(length(sample.names))){
  sample <- sample.names[i]
  svs <- all.sv.samples[[sample]][all.sv.samples[[sample]]$SVTYPE=="ICR"]
  names(svs) <- paste0(seqnames(svs),":",start(svs),"-",end(svs))
  
  transloc <-c(transloc,svs)
  transloc$SVTYPE <- "ICR"
  
  svs <- all.sv.samples[[sample]][all.sv.samples[[sample]]$SVTYPE=="INV"]
  names(svs) <- paste0(seqnames(svs),":",start(svs),"-",end(svs))
  inv <- c(inv,svs)
  inv$SVTYPE <- "INV"
  
  svs <- all.sv.samples[[sample]][all.sv.samples[[sample]]$SVTYPE=="DUP"]
  names(svs) <- paste0(seqnames(svs),":",start(svs),"-",end(svs))
  dup <- c(dup,svs)
  dup$SVTYPE <- "DUP"
  
  svs <- all.sv.samples[[sample]][all.sv.samples[[sample]]$SVTYPE=="DEL"]
  names(svs) <- paste0(seqnames(svs),":",start(svs),"-",end(svs))
  del <- c(del,svs)
  del$SVTYPE <- "DEL"
}


SV <- c(transloc,inv,dup,del)
#Filtering common SVs

SV <- SV[SV$IMPRECISE ==FALSE]
dup.ends <- table(SV$END[SV$SVTYPE %in% c("ICR","INV")][duplicated(SV$END[SV$SVTYPE %in% c("ICR","INV")])])
todel <- SV$END[SV$SVTYPE %in% c("ICR","INV")][SV$END[SV$SVTYPE %in% c("ICR","INV")]%in%names(dup.ends[dup.ends>=3])]
SV <- SV[!SV$END %in%todel]

dup.ends <- table(names(SV)[duplicated(names(SV))])
todel <- names(SV)[names(SV)%in%names(dup.ends[dup.ends>=3])]
SV <- SV[!SV$END %in%todel]

table(SV$Sample)

#Filtering common SVs

dup.ends <- table(transloc$END[duplicated(transloc$END)])
todel <- transloc[transloc$END%in%names(dup.ends[dup.ends>=3])]
SV <- SV[!SV$END %in%todel]

# table(todel$END)
# tofilter <- transloc[transloc$END%in%names(dup.ends[dup.ends<3])]
# table(tofilter$END)

transloc <- transloc[!names(transloc)%in%names(todel)]
dup.starts <- table(names(transloc)[duplicated(names(transloc))])
todel <- transloc[names(transloc)%in%names(dup.starts[dup.ends>=3])]
transloc <- transloc[!names(transloc)%in%names(todel)]

table(transloc$Sample)
# length(tt)
# length(transloc)
# table(tt$END[duplicated(tt$END)])
# todel[todel$END =="chr8:51817583-51817583"]
#
# transloc <- transloc[transloc$END%in%names(table(transloc$END)[table(transloc$END)<3])]
# transloc <- transloc[names(table(transloc$END)[table(transloc$END)<3])]


dup.ends <- table(inv$END[duplicated(inv$END)])
todel <- inv[inv$END%in%names(dup.ends[dup.ends>=3])]
inv <- inv[!names(inv)%in%names(todel)]
dup.start <- table(names(inv)[duplicated(names(inv))])

todel <- inv[names(inv)%in%names(dup.start[dup.start>=3])]
inv <- inv[!names(inv)%in%names(todel)]

# inv <- inv[inv$END%in%names(table(inv$END)[table(inv$END)<3])]

# inv <- inv[!duplicated(inv$END)]

#DUP
# dup.ends <- table(dup$END[duplicated(dup$END)])
# todel <- dup[dup$END%in%names(dup.ends[dup.ends>=3])]
# dup <- dup[!names(dup)%in%names(todel)]
dup.start <- table(names(dup)[duplicated(names(dup))])
todel <- dup[names(dup)%in%names(dup.start[dup.start>=3])]
dup <- dup[!names(dup)%in%names(todel)]

#DEL
dup.start <- table(names(del)[duplicated(names(del))])
todel <- del[names(del)%in%names(dup.start[dup.start>=3])]
del <- del[!names(del)%in%names(todel)]
table(dup.start)

table(SV$Sample)

SV <- c(transloc,inv,dup,del)

# df.fusion <- read.table("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_Tumors/WGS/results/TCGA_ChiTaRS_combined_fusion_information_on_hg19.txt",sep = "\t", header = F )
# 
# gr.bpoint1 <- toGRanges(paste0(df.fusion$V6,":",df.fusion$V7,"-",df.fusion$V7))
# names(gr.bpoint1)<- df.fusion$V4
# gr.bpoint2 <- toGRanges(paste0(df.fusion$V10,":",df.fusion$V11,"-",df.fusion$V11))
# names(gr.bpoint2)<- df.fusion$V4
# 
# chain <- import.chain("/imppc/labs/eslab/mmagallon/Annotatiions/Liftover_chains/hg19ToHg38.over.chain")
# liftover.bp1 <- unlist(liftOver(x = gr.bpoint1, chain = chain))
# liftover.bp2 <- unlist(liftOver(x = gr.bpoint2, chain = chain))
# 
# ss <- subsetByOverlaps(transloc,liftover.bp2,type = "any")
# ss2 <- subsetByOverlaps(liftover.bp2,ss)
# df.fusion[df.fusion$V4 %in% unique(names(ss2)),]
# 
# ss <- subsetByOverlaps(transloc,liftover.bp2)
# ss2 <- subsetByOverlaps(liftover.bp2,ss)
# df.fusion[df.fusion$V4 %in% unique(names(ss2)),]

#Filttering common CNVs annotated in UCSC
for (s in seq_len(length(sample.names))){
  sn <- sample.names[s]
  sv <- SV[SV$Sample == sn]
  
  #find overlaps to annotate
  mcols(sv)$DGV_symbol <-""
  mcols(sv)$DGV_SV <-""
  str.end <- c("start", "end")
  
  for(i in seq_len(length(str.end))){
    fnd <- findOverlaps(query = cnv.ucsc.hg38, sv, type = str.end[i])
    fnd <- data.frame(fnd)
    
    sv[fnd$subjectHits]$DGV_symbol <- cnv.ucsc.hg38[fnd$queryHits]$V4
    sv[fnd$subjectHits]$DGV_SV <- cnv.ucsc.hg38[fnd$queryHits]$V9
  }
  
  #repeatMasker
  sv$repeatMasker <- ""
  
  fnd <- findOverlaps(query = toGRanges(rptmskr),sv,type = "any")
  fnd <- data.frame(fnd)
  sv$repeatMasker[fnd$subjectHits] <- rptmskr$V4[fnd$queryHits]
  all.sv.samples[[sn]] <- sv
  
}


#Filtering of common SVs
svs <- unlist(as(all.sv.samples, "GRangesList"))
names(svs) <- paste0(as.character(seqnames(svs)), ":", start(svs),"-", end(svs))
#duplicated SVs
dup.svs <- names(svs)[which(duplicated(svs))]

##### Filtering and gene annotation of the variants
i=1
filt.list <- list()
for(i in seq_len(length(all.sv.samples))){
  sample <- names(all.sv.samples)[i]
  lumpy.dir <- paste0(execution.dir,"/results/",sample,"/Lumpy")
  
  message(sample, " Lumpy Filtering...")
  sv.data <- all.sv.samples[[sample]]
  #Deletion of common SVs
  sv.data <- sv.data[!names(sv.data) %in% dup.svs]
  sv.data <- toDataframe(sv.data)
  
  #Annotation of SV
  loc <- VariantAnnotation::locateVariants(query = toGRanges(sv.data),txdb,AllVariants())
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
  
  sv.data$Genes <- ""
  for(l in seq_len(length(unique(names(loc))))){
    n.loc<- unique(names(loc))[l]
    sv.data[rownames(sv.data) == n.loc,"Genes"] <- paste0(unique(loc$SYMBOL[names(loc) == n.loc]),collapse =",")
    
  }
  
  #Annotation of mate SV ICR
  icr <- toGRanges(sv.data[sv.data$SVTYPE == "ICR","END"])
  names(icr) <- rownames(sv.data[sv.data$SVTYPE == "ICR",])
  
  loc <- VariantAnnotation::locateVariants(query =  icr, txdb, AllVariants())
  
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
  
  sv.data$Mate.Genes <- ""
  for(l in seq_len(length(unique(names(loc))))){
    n.loc<- unique(names(loc))[l]
    sv.data[rownames(sv.data) == n.loc,"Mate.Genes"] <- paste0(unique(loc$SYMBOL[names(loc) == n.loc]),collapse =",")
    
  }
  sv.data$Position <- rownames(sv.data)
  
  write.table(sv.data,file = file.path(lumpy.dir,paste0(sample,"_Lumpy_AllRegions_CommonSVinSamples.csv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  sv <- toGRanges(sv.data, genome)
  
  filt.list[[sample]][["AllRegions"]] <- sv
  
  sv.data.filt <- sv.data[sv.data$DGV_SV == "",]
  write.table(sv.data.filt,file = file.path(lumpy.dir,paste0(sample,"_Lumpy_AllRegions_FilteredCommonSVinSamples.csv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  sv.data.filt <- toGRanges(sv.data.filt, genome)
  
  filt.list[[sample]][["AllRegionsFiltered"]] <- sv.data.filt
  sv.sel <- subsetByOverlaps(sv.data.filt,gene.markers.gf)
  sv.sel.data <- toDataframe(sv.sel)
  sv.sel.data$Position <- rownames(sv.sel.data)
  write.table(sv.sel.data, file.path(lumpy.dir,paste0(sample,"_Lumpy_SV_Selected_Regions_FilteredCommonSVinSamples.csv")),sep = "\t", col.names = TRUE, row.names = FALSE)
  filt.list[[sample]][["SelectedRegionsFiltered"]] <- sv.sel
  
}
save(filt.list,file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.2.RData")))
# save(filt.list,file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.RData")))
# load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.RData")))
i =1
######### Fusion genes DDBB ##### 
filt.list$`HS-Sch-2`
fus <- subsetByOverlaps(all.fusion.genes.cosmic.gr, filt.list$`HS-Sch-2`$AllRegionsFiltered[filt.list$`HS-Sch-2`$AllRegionsFiltered$SVTYPE =="INV"])
ss <-data.frame(fus)

#List of fusion genes from cosmic gene list
fusion.gene <- sort(all.genes.cosmic.gr[grepl("fusion",all.genes.cosmic.gr$Role.in.Cancer)])
i=8
for(i in seq_len(length(sample.names))){
  sn <- sample.names[i]
  gr.data <- filt.list[[sn]]$AllRegionsFiltered
  fus.genes <- subsetByOverlaps(gr.data, fusion.gene )

  bkpnt.fus <- fus.genes[fus.genes$SVTYPE %in% c("INV","ICR")]
  data.frame(bkpnt.fus)
  fus.geness <- subsetByOverlaps(fusion.gene, bkpnt.fus)
  bkpnt.fus <-bkpnt.fus[bkpnt.fus$Genes %in% fus.geness$Gene.Symbol]
  fus.geness <-fus.geness[fus.geness$Gene.Symbol%in%bkpnt.fus$Genes ]
  
}

fus.genes <- subsetByOverlaps(filt.list$`HS-Sch-2`$AllRegionsFiltered[filt.list$`HS-Sch-2`$AllRegionsFiltered$SVTYPE =="ICR"],all.fusion.genes.cosmic.gr )
fus.genes$Position

fusion.gene <- sort(all.genes.cosmic.gr[grepl("fusion",all.genes.cosmic.gr$Role.in.Cancer)])
fus.genes <- subsetByOverlaps(filt.list$`HS-Sch-2`$AllRegionsFiltered,fusion.gene )
fus.genes[fus.genes$SVTYPE %in% c("INV")]
data.frame(fus.genes[fus.genes$SVTYPE %in% c("INV")])
fus.geness <- subsetByOverlaps(fusion.gene ,filt.list$`HS-Sch-2`$AllRegionsFiltered[filt.list$`HS-Sch-2`$AllRegionsFiltered$SVTYPE %in% c("INV")])


df.fusion <-read.table("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS/TCGA_ChiTaRS_combined_fusion_information_on_hg19.txt", sep = "\t", header = F)
gr.bpoint1 <- toGRanges(paste0(df.fusion$V6,":",df.fusion$V7,"-",df.fusion$V7))
names(gr.bpoint1)<- df.fusion$V4
gr.bpoint2 <- toGRanges(paste0(df.fusion$V10,":",df.fusion$V11,"-",df.fusion$V11))
names(gr.bpoint2)<- df.fusion$V4

chain <- import.chain("/imppc/labs/eslab/mmagallon/Annotatiions/Liftover_chains/hg19ToHg38.over.chain")
liftover.bp1 <- unlist(liftOver(x = gr.bpoint1, chain = chain))
liftover.bp2 <- unlist(liftOver(x = gr.bpoint2, chain = chain))

ss <- subsetByOverlaps(filt.list$`HS-Sch-2`$AllRegionsFiltered,liftover.bp2,type = "any")
ss[ss$SVTYPE=="INV"]
ss2 <- subsetByOverlaps(liftover.bp2,ss[ss$SVTYPE=="ICR"])
df.fusion[df.fusion$V4 %in% unique(names(ss2)),]

ss <- subsetByOverlaps(transloc,liftover.bp2)
ss2 <- subsetByOverlaps(liftover.bp2,ss)
df.fusion[df.fusion$V4 %in% unique(names(ss2)),]


filt.list$`HS-Sch-2`$AllRegionsFiltered
