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
library(BSgenome.Hsapiens.UCSC.hg38)
###################### Functions ############################
source(file = "./LumpySVsFunctions.R")

###################### Parameters ############################

#### Cell Lines ####
# Directories
execution.dir <- "./MPNST_cellLines/WGS"
n <- 2 # Number of the same breakpoint to delete

# Sample information
sample.names <- read.table(file.path(execution.dir,"Sample.info.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sample.names <- sample.names$Sample.name

#  Lumpy file
vcf.file <-file.path(execution.dir,"results/Lumpy_smoove_hg38/MPNST_cellLines_mw6-smoove.genotyped.vcf.gz")

# Lumpy dir
for(sample in sample.names){
  # images.dir <- file.path(execution.dir,"results", sample, "CNVkit/ExcludedRegions2","images")
  lumpy.dir <- file.path(execution.dir,"results", sample, "Lumpy")
  if(!file.exists(lumpy.dir))dir.create(lumpy.dir)
  # if(!file.exists(images.dir))dir.create(images.dir)
}

# Lumpy file to save
lumpy.file <- "All_Lumpy_SV_Regions_FilteredCommonSVinSamples.v4.RData"

#### Tumors ####
# # Directories
# execution.dir <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS"
# n <- 4 # Number of the same breakpoint to delete
#
# # Sample information
# sample.names <- read.table(file = file.path(execution.dir,"sample.info.tumors.csv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# # sample.names <- sample.names$Sample.name[c(2,3,8,12)]
# sample.names <- sample.names$Sample.name[c(6,13:nrow(sample.names))]
#
# # lumpy file
# vcf.file <-file.path(execution.dir,"results/LumpySmoove_TumorNormal/SJdD_TumorsNormal_mw10-smoove.genotyped.vcf.gz")
# # Lumpy dir
# for(sample in sample.names){
#   # images.dir <- file.path(execution.dir,"results", sample, "CNVkit/ExcludedRegions2","images")
#   lumpy.dir <- file.path(execution.dir,"results", sample, "Lumpy")
#   if(!file.exists(lumpy.dir))dir.create(lumpy.dir)
#   # if(!file.exists(images.dir))dir.create(images.dir)
# }
# 
# # Lumpy file to save
# lumpy.file <- "All_Lumpy_SV_Regions_FilteredCommonSV_SJdD_Samples.v4.RData"


###### Loading iformation from data bases needed to filter the variants #####

# Genome
ref.genome.fq <- "/imppc/labs/eslab/mmagallon/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
ref.genome <- "hg38"
orgdb <- org.Hs.eg.db
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- "hg38"

# repeatMasker 
rptmskr <- read.table(file.path("/imppc/labs/eslab/mmagallon/Annotatiions/RepeatMasker/rmskr_hg38.bed"),header = FALSE, stringsAsFactors = FALSE)
rptmskr <- toGRanges(rptmskr, genome = ref.genome)
rptmskr <-filterChromosomes(rptmskr)

# CNV UCSC track
cnv.ucsc.hg38 <- read.table("/imppc/labs/eslab/mmagallon/Annotatiions/DGV_population.bed", sep = "\t",header = FALSE)
cnv.ucsc.hg38 <- toGRanges(cnv.ucsc.hg38)
unique(cnv.ucsc.hg38$V9)
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "139,69,19"] <- "Loss_Gain"
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "200,0,0"] <- "Loss"
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "200,0,200"] <- "Inv"
cnv.ucsc.hg38$V9[cnv.ucsc.hg38$V9 == "0,0,200"] <- "Gain"
cnv.ucsc.hg38 <- sortSeqlevels(cnv.ucsc.hg38)
cnv.ucsc.hg38 <- sort(cnv.ucsc.hg38)

# Selection of interesting genes
gene.markers <- read.table(file.path("./special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
gene.markers.gf <- toGRanges(gene.markers$V2,genome = genome)
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6
regions.df <- toDataframe(regions.selected)
gene.markers.gf <- gene.markers.gf[gene.markers.gf$Genes !="BCR"]

# TSG  COSMIC list: this was downloaded from COSMIC website but make reference to uk.sanger: http://cancer.sanger.ac.uk/census,
# this is de webpage to download this information 
# https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v95/cancer_gene_census.csv?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1642007506&Signature=oVGW10M4sf3rGTeU06M1bSFFlys%3D

all.genes.cosmic <- read.table(file.path("./Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]

tumor.sup.gr <- all.genes.cosmic.gr[grepl("TSG", all.genes.cosmic.gr$Role.in.Cancer)]
# tumor.sup.gr <- toGRanges(tumor.sup$Genome.Location)
# mcols(tumor.sup.gr) <- tumor.sup[,c(1,2,5,6:20)]

#### Fusion Genes Cosmic ####
all.fusion.genes.cosmic <- read.csv(file.path("./CosmicFusionExport.tsv"), sep = "\t",header = T,comment.char = "")

all.fusion.genes.cosmic <- data.frame(chr = paste0("chr",all.fusion.genes.cosmic$X3._CHROMOSOME),
                                      start = all.fusion.genes.cosmic$X3._GENOME_START_FROM,
                                      end = all.fusion.genes.cosmic$X3._GENOME_STOP_FROM,
                                      gene = all.fusion.genes.cosmic$X3._GENE_NAME,
                                      fusion.gene = all.fusion.genes.cosmic$X5._GENE_NAME,
                                      fusion = all.fusion.genes.cosmic$TRANSLOCATION_NAME,
                                      histology = all.fusion.genes.cosmic$PRIMARY_HISTOLOGY)
all.fusion.genes.cosmic <- all.fusion.genes.cosmic[all.fusion.genes.cosmic$chr !="chrNA",]
all.fusion.genes.cosmic.gr <- toGRanges(all.fusion.genes.cosmic, genome = "hg38")




############### Loading Lumpy results ############

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
  all.sv.samples[[sample.name]] <- filterChromosomes(sort(svs))
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
table(SV$Sample)
#Filtering common SVs
# todel <- names(table(names(SV)[duplicated(SV$END)])[table(names(SV)[duplicated(SV$END)])>=2])
# SV <- SV[!SV$END %in%todel]
# todel <- names(table(names(SV)[duplicated(names(SV))])[table(names(SV)[duplicated(names(SV))])>=2])
# SV <- SV[!names(SV) %in% todel]
# length(SV)
# SV <- SV[SV$IMPRECISE ==FALSE]

# Number of duplicated SV between samples

end <- SV$END
end <- end[duplicated(end)]

dup.rows <- SV[SV$END %in% end]
todel.1 <- table(dup.rows$END)>=n
todel.1 <- names(todel.1[todel.1 == TRUE])
todel.1 <- names(SV)[SV$END %in% todel.1]

start <- names(SV)
start <- start[duplicated(start)]
dup.rows <- SV[names(SV) %in% start]
todel.2 <- table(names(dup.rows))>=n
todel.2 <- names(todel.2[todel.2==TRUE])
todel <- unique(c(todel.1,todel.2))
length(todel)
length(SV)

SV.filt <- SV[!names(SV) %in% todel]
length(SV.filt)

SV.filt <- SV.filt[SV.filt$IMPRECISE ==FALSE]

# Filttering common CNVs annotated in UCSC
SV <- c()
for (s in seq_len(length(sample.names))){
  sn <- sample.names[s]
  sv <- SV.filt[SV.filt$Sample == sn]
  
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
  SV <- c(SV,sv)
}
SV <- unlist(as(SV, "GRangesList"))


##### Filtering variants annotated
filt.list <- list()
for (i in seq_len(length(sample.names))){
  sample <- sample.names[i]
  lumpy.dir <- paste0(execution.dir,"/results/",sample,"/Lumpy")
  if(!file.exists(lumpy.dir)) dir.create(lumpy.dir)
  sv.data <- toDataframe(unique(SV[SV$Sample == sample]))
  message("Annotating SVs from ", sample)
  #Annotation of SV
  sv <- SVsAnnotation(sv.data = sv.data)
  
  #Saving results without filtering
  sv.data <- toDataframe(sv)
  # write.table(sv.data,file = file.path(lumpy.dir,paste0(sample,"_Lumpy_AllRegions_CommonSVinSamples.csv")), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  filt.list[[sample]][["AllRegions"]] <- sv
  
  sv.data.filt <- sv.data[sv.data$DGV_SV == "",]
  # Saving data filtered
  # write.table(sv.data.filt,file = file.path(lumpy.dir,paste0(sample,"_Lumpy_AllRegions_FilteredCommonSVinSamples.csv")), sep = "\t", col.names = TRUE, row.names = FALSE)
   
  sv.data.filt <- toGRanges(sv.data.filt, genome)
  
  filt.list[[sample]][["AllRegionsFiltered"]] <- sv.data.filt
  
  # Selecting SVs from specific genes 
  sv.sel <- subsetByOverlaps(sv.data.filt, gene.markers.gf)
  sv.sel.data <- toDataframe(sv.sel)
  sv.sel.data$Position <- rownames(sv.sel.data)
  # Saving data of specific genes
  # write.table(sv.sel.data, file.path(lumpy.dir,paste0(sample,"_Lumpy_SV_Selected_Regions_FilteredCommonSVinSamples.csv")),sep = "\t", col.names = TRUE, row.names = FALSE)
  
  filt.list[[sample]][["SelectedRegionsFiltered"]] <- sv.sel

}

# save(filt.list, file = file.path(execution.dir,"results/", lumpy.file))
load(file = file.path(execution.dir,"results/", lumpy.file))
i =1

# ######### Fusion genes DDBB ##### 
# filt.list$`HS-Sch-2`
# fus <- subsetByOverlaps(all.fusion.genes.cosmic.gr, filt.list$`HS-Sch-2`$AllRegionsFiltered[filt.list$`HS-Sch-2`$AllRegionsFiltered$SVTYPE =="INV"])
# ss <-data.frame(fus)
# 
# #List of fusion genes from cosmic gene list
# fusion.gene <- sort(all.genes.cosmic.gr[grepl("fusion",all.genes.cosmic.gr$Role.in.Cancer)])
# i=4
# for(i in seq_len(length(sample.names))){
#   sn <- sample.names[i]
#   gr.data <- filt.list[[sn]]$AllRegions
#   if(is.null(gr.data)) next
#   gr.data <- gr.data[gr.data$SVTYPE %in% c("INV","ICR")]
#   
#   bkpnt.fus <- subsetByOverlaps(gr.data, fusion.gene)
#   # bkpnt.fus.2 <- subsetByOverlaps(toGRanges(gr.data$END,genome=genome),fusion.gene)
#  
#   
#   data.frame(bkpnt.fus)
#   fus.geness <- subsetByOverlaps(fusion.gene, bkpnt.fus)
#   # fus.geness <- subsetByOverlaps(fusion.gene, bkpnt.fus.2)
#   
#   bkpnt.fus <-bkpnt.fus[bkpnt.fus$Genes %in% fus.geness$Gene.Symbol]
#   fus.geness <-fus.geness[fus.geness$Gene.Symbol%in%bkpnt.fus$Genes ]
#   print(bkpnt.fus)
#   print(fus.geness)
# }
# 
