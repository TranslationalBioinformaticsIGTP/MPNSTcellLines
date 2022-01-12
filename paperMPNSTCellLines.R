#############################################################################
#         Script to process data to generate figures of the papaer          #
#############################################################################
# Packages Needed
library(CopyNumberPlots)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)

######################## Functions #########################################
source(file = "/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/loadingLumpySVs.R")
source(file = "./MPNST_cellLines/funtionsPaperCellLines.R")
################## Parameters #######################
execution.dir <- "./MPNST_cellLines/WGS"

#Cosmic and Gene files
#Selection of interesting genes
gene.markers <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_cellLines/WGS","special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
gene.markers.gf <- toGRanges(gene.markers$V2,genome = genome)
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6
regions.df <- toDataframe(regions.selected)
gene.markers.gf <- gene.markers.gf[gene.markers.gf$Genes !="BCR"]
gene.markers.gf<- gene.markers.gf[gene.markers.gf$Genes %in% c("NF1","SUZ12","EED",
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
                                                       # "NF2",
                                                       # "LZTR1",
                                                       "CXCR4",
                                                       "TWIST1"
)]
# TSG  COSMIC list
all.genes.cosmic <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS","Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]

all.genes.cosmic[545,]
tumor.sup.gr <- all.genes.cosmic.gr[grepl("TSG", all.genes.cosmic.gr$Role.in.Cancer)]



#Fusion Genes Cosmic

table(all.fusion.genes.cosmic$PRIMARY_HISTOLOGY)
all.fusion.genes.cosmic <- read.csv(file.path("./MPNST_tumors/WGS","CosmicFusionExport.tsv"), sep = "\t",header = T,comment.char = "")

all.fusion.genes.cosmic <- data.frame(chr = paste0("chr",all.fusion.genes.cosmic$X3._CHROMOSOME),
                                      start = all.fusion.genes.cosmic$X3._GENOME_START_FROM,
                                      end = all.fusion.genes.cosmic$X3._GENOME_STOP_FROM,
                                      gene = all.fusion.genes.cosmic$X3._GENE_NAME,
                                      fusion.gene = all.fusion.genes.cosmic$X5._GENE_NAME,
                                      fusion = all.fusion.genes.cosmic$TRANSLOCATION_NAME,
                                      histology = all.fusion.genes.cosmic$PRIMARY_HISTOLOGY)
all.fusion.genes.cosmic <- all.fusion.genes.cosmic[all.fusion.genes.cosmic$chr !="chrNA",]
all.fusion.genes.cosmic.gr <- toGRanges(all.fusion.genes.cosmic, genome = "hg38")


# Loading Cell lines's sample info
sample.data <- read.table(file.path(execution.dir,"Sample.info.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)[-c(9,10),]
sample.names <- sample.data$Sample.name
sample.names <- sample.names[c(3,2,5,1,4,6,8,7)]
name.cells <- c("90-8-TL","ST88-14", "NMS-2", "S462", "sNF96.2",  "STS-26T", "HS-PSS", "HS-Sch-2" )

names(sample.names) <- name.cells
sample.names <- sample.names[c(4,2,1,5,3,6,8,7)]

#other parameters
orgDb <- org.Hs.eg.db
genome <- "hg38"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Directories
for(sample in sample.names){
  # images.dir <- file.path(execution.dir,"results", sample, "CNVkit/ExcludedRegions2_ploidy_Thr","images")
  images.dir <- file.path(execution.dir,"results",sample,"images")
  lumpy.dir <- file.path(execution.dir,"results", sample, "Lumpy")
  if(!file.exists(lumpy.dir))dir.create(lumpy.dir)
  if(!file.exists(images.dir))dir.create(images.dir)
}

############ Loading and filtering of the data to plot ####################
# Loading the results of CNVkit
sample.log <- list()
# sample.seg <- list()
sample.cn <- list()
sample <- "S462"
i=4
for(i in seq_len(length(sample.names))){
  sample <- sample.names[i]
  sn <- names(sample.names)[i]
  # Data to plot
  # logs <-loadCopyNumberCallsCNVkit(file.path(execution.dir,"results", sample,"CNVkit/ExcludedRegions2", paste0(sample,".cnr")))
  logs <-loadCopyNumberCallsCNVkit(file.path(execution.dir,"results", sample,"CNVkit/ExcludedRegions2_ploidy", paste0(sample,".cnr")))
  
  names(mcols(logs))[2] <- "lrr"
  sample.log[[sn]] <-logs
  # sample.seg[[sample]] <- loadCopyNumberCallsCNVkit(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_Tumors/WGS/results",sample,"CNVkit/ExcludedRegions", paste0(sample,"_seg.cns")))
  # sample.cn[[sample]] <- loadCopyNumberCallsCNVkit(file.path(execution.dir, "results", sample, "CNVkit/ExcludedRegions2", paste0(sample,"_call_4.2.cns")))
  sample.cn[[sn]] <- loadCopyNumberCallsCNVkit(file.path(execution.dir, "results", sample, "CNVkit/ExcludedRegions2_ploidy", paste0(sample,"_call_ploidy.cns")))
  
  # loh
  loh <- sample.cn[[sn]]$cn1 != sample.cn[[sn]]$cn2
  sample.cn[[sn]]$loh <-loh
}
# save(sample.log, file= file.path(execution.dir, "results/CNVkit_lrr_AllSamples.RData"))
# save(sample.cn, file= file.path(execution.dir, "results/CNVkit_cn_AllSamples.RData"))


#loading and filteringcliffHunteR Data
cliff.info.samples <- c()
for(i in seq_len(length(sample.names))){
  sample <- sample.names[i]
  result.dir <- file.path(execution.dir,"results",sample)
  cliff.dir <- file.path(result.dir,"CliffHunteR")
  cliff.data <- read.table(file.path(cliff.dir, paste0(sample,"_CliffHunter_CNV_0.25_4.csv")),sep = "\t", header = TRUE,check.names = T)
  cliff.data <- toGRanges(cliff.data,genome = genome)
  cliff.data <- cliff.data[cliff.data$percentage >=0.6]
  cliff.data$sample <- sample
  cliff.info.samples <- c(cliff.info.samples,cliff.data)

}

cliff.info.samples <-unlist(as(cliff.info.samples, "GRangesList"))
names(cliff.info.samples) <- paste0(seqnames(cliff.info.samples),":",start(cliff.info.samples),"-", end(cliff.info.samples))
#Filtering duplicated SVs. Artefacts
cliff.info.samples <-cliff.info.samples[!duplicated(names(cliff.info.samples))]
#Split a grnages by sample
cliff.info.samples <- split(cliff.info.samples, cliff.info.samples$sample)
# cliff.info.samples$S462
# cliff.info.s462 <- subsetByOverlaps(cliff.info.samples$S462, gene.markers,type = "within")
cliff.info.sample.genes <- lapply(cliff.info.samples,function(x){subsetByOverlaps(x, gene.markers.gf)})

#Lumpy Filtered Data
# load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.RData")))
# load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.2.RData")))
# load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.v3.RData")))
load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.v4.RData")))


#Selected genes
lumpy.selgenes  <- lapply(filt.list,function(x){lapply(x,function(y){subsetByOverlaps(y, gene.markers.gf)})})
lumpy.selgenes <- lapply(lumpy.selgenes, function(x){x[["AllRegionsFiltered"]]} )

#Joining CliffHunter and Lumpy results for plotting
join.cliff.lumpy <- list()
i=1
for(i in seq_len(length(sample.names))){
  sn <- sample.names[i]
  lump <- lumpy.selgenes[[sn]]
  lump <-lump [lump$IMPRECISE == FALSE]
  mcols(lump) <- mcols(lump)[c("END","SVTYPE","Genes","Sample")]
  cliff <- c()
  if(length(cliff.info.sample.genes[[sn]])>0){
    cliff <- cliff.info.sample.genes[[sn]]
    mcols(cliff) <- mcols(cliff)[c("mate.pos","SYMBOL","sample")]
    mcols(cliff)$SVTYPE <- ""
    colnames(mcols(cliff)) <- c("END","Genes", "Sample","SVTYPE")
    mcols(cliff) <- mcols(cliff)[c("END","SVTYPE","Genes","Sample")]
  }


  join.cliff.lumpy[[sn]] <- sort(unique(c(lump,cliff)))

}
#After joining the data, we manually filtered the false positive variants comming from cliff data.

join.cliff.lump.sn <- sort(unlist(as(join.cliff.lumpy,"GRangesList")))
names(join.cliff.lump.sn)<-gsub(pattern = "NMS.2",replacement = "NMS-2", names(join.cliff.lump.sn))
names(join.cliff.lump.sn)<-gsub(pattern = "SNF96.2",replacement = "SNF96-2", names(join.cliff.lump.sn))

nms <- strsplit(x = names(join.cliff.lump.sn),"\\.")
sn <- unlist(lapply(nms, function(x){x[1]}))
nms <- unlist(lapply(nms, function(x){x[2]}))
names(join.cliff.lump.sn) <- nms


# join.cliff.lumpy$`90-8TL` <- join.cliff.lumpy$`90-8TL`[-c(2,3,5,6,7,8,9)]

join.cliff.lumpy$`90-8TL`$SVTYPE
toDataframe(join.cliff.lumpy$`ST88-14`)

join.cliff.lumpy$`ST88-14`$SVTYPE[11] <- "ICR"
join.cliff.lumpy$`ST88-14` <- join.cliff.lumpy$`ST88-14`[-c(1,2,7,8,9,10)]

toDataframe(join.cliff.lumpy$NMS.2)
join.cliff.lumpy$NMS.2$SVTYPE[2] <- "ICR"
join.cliff.lumpy$NMS.2$SVTYPE[8] <- "ICR"
join.cliff.lumpy$NMS.2 <- join.cliff.lumpy$NMS.2[join.cliff.lumpy$NMS.2$SVTYPE%in%c("ICR","INV")]

toDataframe(join.cliff.lumpy$S462)

join.cliff.lumpy$S462 <- join.cliff.lumpy$S462[-c(1,2,3,6:12)]
join.cliff.lumpy$S462$SVTYPE <- c("ICR","ICR","DEL")
names(join.cliff.lumpy$S462)[3] <- "chr17:31940127-31940427"
join.cliff.lumpy$S462$END[3] <- 31940427
end(join.cliff.lumpy$S462)[3] <- 31940427
start(join.cliff.lumpy$S462)[3] <- 31940127



join.cliff.lumpy$`STS-26T` <- join.cliff.lumpy$`STS-26T`[-c(1,2)]

join.cliff.lumpy$`HS-PSS` <- join.cliff.lumpy$`HS-PSS`[-c(1:6,8,9)]

join.cliff.lumpy$`HS-Sch-2` <- join.cliff.lumpy$`HS-Sch-2`[-c(1,2,5)]

# save(join.cliff.lumpy, file= file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.RData"))
# save(join.cliff.lumpy, file= file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.v4.RData"))

# sv.data <- toDataframe(unlist(as(join.cliff.lumpy,"GRangesList")))
# write.table(sv.data, file = file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.v4.csv"),sep = "\t", col.names = TRUE,row.names = FALSE)


#Strelka Filtered Variants
# load(file = file.path(execution.dir, "results/All_samples_Strelka_Variants.RData"))
# load(file = file.path(execution.dir, "results/All_samples_Strelka_Variants_NewFilters.RData"))
load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.3.RData"))

strelka.files <- strelka.variants
# lapply(strelka.files,length)
# strelka.files
strk <- unlist(as(strelka.files,"GRangesList"))
for(i in seq_len(length(sample.names))){
  names(strk)[grepl(sample.names[i],names(strk))] <- gsub(paste0(sample.names[i],"."),"",names(strk)[grepl(sample.names[i],names(strk))])
}  
strk$Change <- names(strk)

todel <- table(names(strk))[table(names(strk))>=2]
strk <- strk[!names(strk)%in%names(todel)]


# Cosmic genes
strelka.var <-subsetByOverlaps(strk, all.genes.cosmic.gr)

for(i in seq_len(length(sample.names))){
  sn <- sample.names[i]
  strelka.sn.cos <- strelka.var[strelka.var$Sample ==sn]
  strelka.df <- toDataframe(sort(strelka.sn.cos))
  # write.table(strelka.df, file.path(execution.dir,"results",sn,"Strelka_GermVariants2.9.10/results", paste0(sn,"_SNVs_cosmicGenes.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(strelka.df, file.path(execution.dir,"results",sn,"Strelka_GermVariants2.9.10/results", paste0(sn,"_SNVs_cosmicGenes.toKeep.3.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}

# Selected genes 
strelka.var <-subsetByOverlaps(strk, gene.markers.gf)

strelka.var[unlist(strelka.var$Func.refGene) =="splicing"]
mcols(strelka.var) <- mcols(strelka.var)[c("Gene.refGene",
                                           "ExonicFunc.refGene",
                                           "Func.refGene",
                                           "AAChange.refGene",
                                           "CLNSIG","SIFT_pred","Polyphen2_HDIV_pred",
                                           "Polyphen2_HVAR_pred",
                                           "LRT_pred","FATHMM_pred","MutationTaster_pred",
                                           "MutationAssessor_pred","InterVar_automated",
                                           "AF","AF_popmax","prop_pathogenic_pred","Ref", "Alt","Sample")]

# strelka.var <- strelka.var[!grepl("Benign",unlist(strelka.var$CLNSIG)) ]
strelka.var$Variant.type <- "SNV"
strelka.var.df <- toDataframe(strelka.var)
strelka.var.df$Changes <- ""
strelka.var.df$Changes <-rownames(strelka.var.df)
strelka.var.df[strelka.var.df$Gene.refGene=="CDKN2A"]
# write.table(strelka.var.df,file.path(execution.dir,"results/SNVs.selectedGenes.AllSamples.csv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(strelka.var.df,file.path(execution.dir,"results/SNVs.selectedGenes.AllSamples.toKeep.3.2.csv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

join.cliff.lumpy <- unlist(as(join.cliff.lumpy,"GRangesList"))
join.cliff.lumpy$Variant.type <-"SV" 

complete.table.points <- c(granges(strelka.var),granges(join.cliff.lumpy))
complete.table.points$Variant <- c(strelka.var$Variant.type, join.cliff.lumpy$Variant.type)
complete.table.points$Variant.type <- c(strelka.var$ExonicFunc.refGene, join.cliff.lumpy$SVTYPE)
complete.table.points$Genes <- c(strelka.var$Gene.refGene,join.cliff.lumpy$Genes)
complete.table.points$Positions <- c(names(strelka.var), join.cliff.lumpy$END)

complete.table.points$Sample <- c(strelka.var$Sample, join.cliff.lumpy$Sample)
complete.table.points.df <- toDataframe(complete.table.points)
write.table(complete.table.points.df,file.path(execution.dir,"results/SNVsandSVs.selectedGenes.AllSamples.toKeep.3.csv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Loading the SNParray data from WGS
load(file.path(execution.dir,"results","MPNST_WGS_cellLines_BAF_Filtered.RData"))

##### loading GAP SNP-array data ####
data.dir <- "/imppc/labs/eslab/bgel/Pipelines/Executions/SNPArrayCopyNumberCalling/v0.2"
result.dir <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/SNP-array/results/paperCDKN2A"

cell.lines <- read.table("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/SNP-array/cellLinesDataBernatDir.csv",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cell.lines <-cell.lines[c(4,2,1,5,3,6,8,7),]
rownames(cell.lines)<- names(sample.names)
sample.names

####Loading segment data
segment.files <- list()
snp.files <- list()
i=4
for (i in seq_len(nrow(cell.lines))){
  #segment data 
  segment.file <- file.path(cell.lines$data.path[i],"results/segments", paste0(cell.lines$result.file[i],".GAP.segments.txt"))
  if(file.exists(segment.file)){
    segment.files[[rownames(cell.lines)[i]]] <- loadCopyNumberCalls(cnv.data = segment.file) 
    segment.files[[rownames(cell.lines)[i]]] <- UCSCStyle(segment.files[[rownames(cell.lines)[i]]])
    segment.files[[rownames(cell.lines)[i]]]$cn[segment.files[[rownames(cell.lines)[i]]]$cn %in% 100] <- 0
    
  } else{
    stop(segment.file, "doesn't exist")
  }
  
  #snpData
  snp.file <- file.path(cell.lines$data.path[i],"data")
  snp.file <- file.path(snp.file, list.files(snp.file))
  snp.files[[rownames(cell.lines)[i]]] <- loadSNPData(snp.data = snp.file) 
  snp.files[[rownames(cell.lines)[i]]] <- UCSCStyle(snp.files[[rownames(cell.lines)[i]]])
}

# Fixing GAP error in LOH
for(i in seq_len(length(segment.files))){
  segment.files[[i]]$loh[segment.files[[i]]$cn ==0] <-0
}


#Gene inactivation table
table.info.plot <- list()
i=4
for(i in seq_len(length(sample.names))){
  sn <- names(sample.names)[i]
  gene.inact <- gene.markers.gf
  complete.list <- complete.table.points[complete.table.points$Sample ==sample.names[i]]
  cn <- subsetByOverlaps(sample.cn[[i]],gene.inact)
  cn$CNV <- ""
  cn$CNV[cn$cn>2] <- "Gain"
  cn$CNV[cn$cn<=1] <- "HetLoss"
  cn$CNV[cn$cn<=0] <- "HomLoss"
  cn$CNV[cn$cn ==2] <- "2n"

  cn.gap <- subsetByOverlaps(segment.files[[i]],gene.inact)
  cn.gap$loh <- as.logical(cn.gap$loh)
  cn.gap$CNV <- ""
  cn.gap$CNV[cn.gap$cn>2] <- "Gain"
  cn.gap$CNV[cn.gap$cn<=1] <- "HetLoss"
  cn.gap$CNV[cn.gap$cn<=0] <- "HomLoss"
  cn.gap$CNV[cn.gap$cn ==2] <- "2n"


  gene.inact$CNV.cnvkit <- ""
  for(j in seq_len(length(cn))){
    pos <- subsetByOverlaps(gene.inact,cn[j])
    if(length(pos)>0){
      gene.inact[names(pos)]$CNV.cnvkit <- cn$CNV[j]

    }
  }


  gene.inact$CNV.GAP <- ""
  gene.inact$LOH.GAP <- ""
  for(j in seq_len(length(cn.gap))){
    pos <- subsetByOverlaps(gene.inact,cn.gap[j])
    if(length(pos)>0){

      gene.inact[names(pos)]$LOH.GAP <- cn.gap$loh[j]
      gene.inact[names(pos)]$CNV.GAP <- cn.gap$CNV[j]
    }
  }

  complete.list$SpecificGene <- ""
  for(j in seq_len(length(gene.inact))){
    pos <- subsetByOverlaps(complete.list,gene.inact[j])
    if(length(pos)>0){
      complete.list[names(pos)]$SpecificGene <- gene.inact$Genes[j]
    }
  }

  gene.inact$Variant <-NA
  gene.inact$Variant.type <- NA
  gene.inact$Start.bkpoint <-NA
  gene.inact$End.bkpoint <- NA
j=17
  for(j in seq_len(length(gene.inact$Genes))){
    gene <- gene.inact$Genes[j]
    variants.inf <- complete.list[complete.list$SpecificGene%in%gene]

    if(length(variants.inf)==0)next

    if(length(variants.inf)==1){
      gene.inact$Variant[j] <- variants.inf$Variant
      gene.inact[j]$Variant.type <- variants.inf$Variant.type
      gene.inact[j]$Start.bkpoint<- paste0(as.character(seqnames(variants.inf)),":",start(variants.inf),"-",end(variants.inf))
      gene.inact[j]$End.bkpoint<- variants.inf$Positions

    }
    if(length(variants.inf)>1){
      rep <- gene.inact[gene.inact$Genes%in%gene]
      for(k in seq_len(length(variants.inf)-1)){
        gene.inact[names(rep)]$Variant<- variants.inf[k]$Variant
        gene.inact[names(rep)]$Variant.type <- variants.inf[k]$Variant.type
        gene.inact[names(rep)]$Start.bkpoint<- paste0(as.character(seqnames(variants.inf[k])),":",start(variants.inf[k]),"-",end(variants.inf[k]))
        gene.inact[names(rep)]$End.bkpoint<- variants.inf$Positions[k]

        rep$Variant <- variants.inf$Variant[k+1]
        rep$Variant.type <- variants.inf$Variant.type[k+1]
        rep$Start.bkpoint<- paste0(as.character(seqnames(variants.inf[k+1])),":",start(variants.inf[k+1]),"-",end(variants.inf[k+1]))
        rep$End.bkpoint<- variants.inf$Positions[k+1]
        gene.inact <- c(gene.inact,rep)
      }
    }


  }

  table.info <- sort(gene.inact)
  names(table.info)<- c(1:length(table.info))

  table.info <-toDataframe(table.info)

  df.del <- table.info[which(table.info$Variant.type == "DEL"),]
  df.del <- df.del[which(df.del$CNV.cnvkit == "Gain"),]

  if(nrow(df.del)>=1){
    thr <- end(complete.list[complete.list$SpecificGene%in%df.del$Genes])- start(complete.list[complete.list$SpecificGene%in%df.del$Genes])
    if(thr >100000){
      table.info <- table.info[!rownames(table.info)%in% rownames(df.del),]

    }else{
      table.info <- table.info
    }

  }

  df.del <- table.info[which(table.info$Variant.type == "DUP"),]
  df.del <- df.del[which(df.del$CNV.cnvkit == "2n"),]

  table.info <- table.info[!rownames(table.info)%in% rownames(df.del),]
  df.del <- table.info[which(table.info$Variant.type == "DEL"),]
  df.del <- df.del[which(df.del$CNV.cnvkit == "2n"),]

  table.info <- table.info[!rownames(table.info)%in% rownames(df.del),]
  df.del <- table.info[which(table.info$Variant.type == "DUP"),]
  df.del <- df.del[which(df.del$CNV.cnvkit == "Gain"),]

  table.info <- table.info[!rownames(table.info)%in% rownames(df.del),]

  table.info$Sample <- sn
  if(sn == "90-8-TL"){
    table.info$Variant.type[table.info$Genes=="TP53"] <- "splicing"
  }
  if(sn == "NMS-2"){
    table.info$Variant.type[table.info$Genes=="NF1"] <-"splicing"

  }
  if(sn == "HS-Sch-2"){
    table.info$Variant.type[table.info$Variant.type == "."] <-"splicing"

  }
  

  # write.table(table.info, file.path(execution.dir,"results",sample.names[i], paste0(sample.names[i],"_SNVsandSVsandCNV.selectedGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(table.info, file.path(execution.dir,"results",sample.names[i], paste0(sample.names[i],"_SNVsandSVsandCNV.selectedGenes.toKeep.3.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  

  table.info.plot[[sn]] <- sort(toGRanges(table.info,genome = genome))

}

df.complete <- toDataframe(sort(unlist(as(table.info.plot,"GRangesList"))))

# write.table(table.info, file.path(execution.dir,"results",sample.names[i], paste0(sample.names[i],"_SNVsandSVsandCNV.selectedGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(df.complete, file.path(execution.dir,"results","ALL_SNVsandSVsandCNV.selectedGenes.toKeep.3.csv"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#TABLE TO PLOT
s = 4
table.to.plot <- list()
for(s in seq_len(length(sample.names))){
  sn <- names(sample.names)[s]
  tp <- table.info.plot[[sn]]

  gene.table <- gene.markers.gf
  gene.table$variant <- ""
  gene.table$variant.type <- ""
  gene.table$cn <- ""
  gene.table$loh <- ""
  gene.table$altered <- ""

  # gene <- "SUZ12"
  i=17
  for(i in seq_len(length(gene.table$Genes))){
    gene <- gene.table$Genes[i]
    if(!gene %in% gene.table$Genes)next
    tp.g <- tp[tp$Genes%in%gene]

    if(!gene %in% tp$Genes){
      tp.g <- c(tp,gene.table[gene.table$Genes == gene])
      tp.g <- tp.g[tp.g$Genes == gene]
    }

    if(length(tp.g)>1){
      variant <- unique(tp.g$Variant)
      variant.type <- unique(tp.g$Variant.type)

      cn <- unique(tp.g$CNV.cnvkit)
      loh <- unique(tp.g$LOH.GAP)

      gene.table[gene.table$Genes == gene]$altered <- "complete"
      gene.table[gene.table$Genes == gene]$loh <- loh
      gene.table[gene.table$Genes == gene]$cn <- cn

      if(length(variant)>1){
        variant <- paste0(variant,collapse = ",")
      }
      if(length(variant.type)>1){
        variant.type <- paste0(variant.type, collapse = ",")
      }
      gene.table[gene.table$Genes == gene]$variant <- variant
      gene.table[gene.table$Genes == gene]$variant.type <- variant.type

    }else{
      gene.table[gene.table$Genes == gene]$variant <- tp.g$Variant
      gene.table[gene.table$Genes == gene]$variant.type <- tp.g$Variant.type
      gene.table[gene.table$Genes == gene]$loh <- tp.g$LOH.GAP
      gene.table[gene.table$Genes == gene]$cn <- tp.g$CNV.cnvkit

      if(tp.g$LOH.GAP ==TRUE & tp.g$Variant%in%c("SV","SNV")){
        gene.table[gene.table$Genes == gene]$altered <- "complete"


      }

      if(tp.g$CNV.cnvkit %in%"HomLoss"){
        gene.table[gene.table$Genes == gene]$altered <- "complete"


      }

      if(tp.g$LOH.GAP ==TRUE & tp.g$CNV.cnvkit %in% "HetLoss"){
        gene.table[gene.table$Genes == gene]$altered <- "complete"


      }
    }

  }
  #Colors gene table
  gene.table$color <- "black"
  gene.table$color[gene.table$altered=="complete"] <- "red"
  gene.table$variant.col <- ""

   #snv
  gene.table$variant.col[grepl("frameshift_insertion", gene.table$variant.type)] <- "grey"
  gene.table$variant.col[grepl("frameshift_deletion", gene.table$variant.type)] <- "black"
  gene.table$variant.col[grepl("frameshift_deletion", gene.table$variant.type)] <- "black"
  gene.table$variant.col[gene.table$variant.type=="splicing"] <- "cornflowerblue"
  gene.table$variant.col[gene.table$variant.type=="nonsynonymous_SNV"] <- "darkgoldenrod3"
  gene.table$variant.col[grepl("stopgain", gene.table$variant.type)] <- "red"

  #sv
  gene.table$variant.col[gene.table$variant.type=="ICR"] <- "black"
  gene.table$variant.col[gene.table$variant.type=="DEL"] <- "green"
  gene.table$variant.col[gene.table$variant.type=="INV"] <- "grey"
  gene.table$variant.col[gene.table$variant.type=="DUP"] <- "brown"



  gene.table$loh.col <- ""
  gene.table$loh.col[gene.table$loh == TRUE] <- 9
  gene.table$Sample <- sn
  table.to.plot[[sn]] <- gene.table
}
table.to.plot$`STS-26T`$altered[table.to.plot$`STS-26T`$Genes=="TP53"]<-"complete"
table.to.plot$`STS-26T`$color[table.to.plot$`STS-26T`$Genes=="TP53"]<-"red"

# data.frame(table.to.plot$S462)
# save(table.info.plot,file = file.path(execution.dir,"results","table.info.plot.2.RData"))
# save(table.to.plot,file = file.path(execution.dir,"results","table.to.plot.toKeep.RData"))
# save(table.info.plot,file = file.path(execution.dir,"results","table.info.plot.RData"))
# save(table.info.plot,file = file.path(execution.dir,"results","table.info.plot.toKeep.RData"))
save(table.to.plot,file = file.path(execution.dir,"results","table.to.plot.toKeep.3.RData"))

# load(file = file.path(execution.dir,"results","table.to.plot.toKeep.3.RData"))
