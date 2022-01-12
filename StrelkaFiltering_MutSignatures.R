###################################################
#          Filtering Strelka variants             #
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
library(VariantAnnotation)

############## Functions ################
source(file = "/imppc/labs/eslab/mmagallon/Projects/cliffhunter/SV.functions.R")

#################### Parameters & Directories #####################
#Directory
execution.dir <-"/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS"
# execution.dir <-"/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS"

#Cell Lines
sample.data <- read.table(file = file.path(execution.dir,"Sample.info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample.names <- sample.data$Sample.name
names(sample.names) <- sample.names

#Tumors
# sample.names <- read.table(file = file.path(execution.dir,"sample.info.tumors.csv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# sample.names <- sample.names$Sample.name[c(2,3,8,12,14,17,22,23)]
# sample.names <- sample.names$Sample.name[c(2,3,8,12)]

ref.genome.fq <- "/imppc/labs/eslab/mmagallon/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref.genome <- "hg38"

#repeatMasker 
# rptmskr <- read.table(file.path("/imppc/labs/eslab/mmagallon/Annotatiions/RepeatMasker/rmskr_hg38.bed"),header = FALSE, stringsAsFactors = FALSE)
# rptmskr <- toGRanges(rptmskr, genome = ref.genome)

gene.markers <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_cellLines/WGS","special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
gene.markers.gf <- toGRanges(gene.markers$V2, genome ="hg38")
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6

# TSG  COSMIC list
all.genes.cosmic <- read.table(file.path("/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS","Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]

# regions.df <- toDataframe(regions.selected)
strelka.files <- list()

################ Filtering parameters ####################
info.param <- c("Gene.refGene",
                "ExonicFunc.refGene",
                "Func.refGene",
                "avsnp150",
                "AAChange.refGene",
                "AF",
                "AF_popmax","GeneDetail.refGene")
geno.param <- "AD"

#Min VAF to filter
min.VAF <- 0.1

for(sn in seq_len(length(sample.names))){
  sample <- sample.names[sn]
  message("Getting Setrelka variants of specific Genes from ", sample)
  # Directories
  result.dir <- file.path(execution.dir,"results",sample)
  strelka.res.path <- file.path(result.dir, "Strelka_GermVariants2.9.10/results/variants", paste0(sample, "_annotated_PASS.vcf.gz.tbi"))


  #Loading strelka variants
  vcf <- VariantAnnotation::readVcf(file = strelka.res.path, genome = ref.genome, param = ScanVcfParam(fixed = NA, geno =geno.param,info=info.param ))
  # vcf <- VariantAnnotation::readVcf(file = strelka.res.path, genome = ref.genome, param = ScanVcfParam(info=info.param ))
  
  
  rowRanges.celltype <- rowRanges(vcf,genome = ref.genome)
  info.tumor <- info(vcf)
  mcols(rowRanges.celltype) <- cbind(mcols(rowRanges.celltype),info.tumor)
  
  #Adding VAF to the sample info
  ad <- VariantAnnotation::geno(vcf)$AD
  ad.ref <- unlist(lapply(ad, "[", 1))
  ad.alt <- unlist(lapply(ad, "[", 2))
  vaf <- ad.alt/(ad.alt + ad.ref)
  vaf[is.nan(vaf)] <- NA
  
  rowRanges.celltype$VAF <- vaf
  
  message(length(rowRanges.celltype)," SNV initials")
  #Unlisting mcols of granges
  for(j in seq_len(ncol(mcols(rowRanges.celltype)))){
    
    if(class(mcols(rowRanges.celltype)[[j]]) =="CompressedCharacterList"){
      ll <- lapply(mcols(rowRanges.celltype)[[j]],function(x)length(x))
      mcols(rowRanges.celltype)[[j]][ll >=2 ]<-lapply(mcols(rowRanges.celltype)[[j]][ll>=2],function(x)paste0(x,collapse = ","))
      mcols(rowRanges.celltype)[j] <- unlist(mcols(rowRanges.celltype)[[j]])
      
    }
  }
  
   message("Starting filtering. Selecting AF< 0.01")
  rowRanges.celltype$AF_popmax <-   as.numeric(rowRanges.celltype$AF_popmax)
  rowRanges.celltype$AF_popmax[is.na(rowRanges.celltype$AF_popmax)] <- "."
  rowRanges.celltypes <- rowRanges.celltype[rowRanges.celltype$AF_popmax<0.01 | rowRanges.celltype$AF_popmax =="."]
  
  message("SNVs after AF filtering ", length(rowRanges.celltypes))
  
  # We delete those variants with UNKNOWN in AAChange.refGene
  rowRanges.celltypes <- rowRanges.celltypes[rowRanges.celltypes$AAChange.refGene !="UNKNOWN"]
  message("SNVs after UNKNOWN filtering ", length(rowRanges.celltypes))
  
  #we also select those varinats with AF = "."
  #If the AF_popmax has a "." we also add this variant as a possible interesting with "."
  # rowRanges.celltypes$Variant_sel_status <- ""
  # rowRanges.celltypes$Variant_sel_status[rowRanges.celltypes$AF_popmax=="."] <- "."
  
  
  seqlevels(rowRanges.celltypes) <- sortSeqlevels(seqlevels(rowRanges.celltypes))
  rowRanges.celltypes <- sort(rowRanges.celltypes)
  rowRanges.celltypes <- filterChromosomes(rowRanges.celltypes)
  strelka.files[[sample]] <-rowRanges.celltypes
}

# save(strelka.files, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignature.3.RData"))
load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignature.3.RData"))
 
# ss <- data.frame(strelka.files$S462)
# 
# After filtering the variants per sample, now we proceed to filter common variants
strk <- unlist(as(strelka.files,"GRangesList"))
names(strk)<-gsub(pattern = "NMS.2",replacement = "NMS-2",names(strk))
names(strk)<-gsub(pattern = "SNF96.2",replacement = "SNF96-2",names(strk))

nms <- strsplit(x = names(strk),"\\.")
sn <- unlist(lapply(nms, function(x){x[1]}))
nms <- unlist(lapply(nms, function(x){x[2]}))
names(strk) <- nms
strk$Sample <- sn
strk$Sample<-gsub(pattern = "SNF96-2",replacement = "SNF96.2",strk$Sample)
strk$Sample<-gsub(pattern = "NMS-2",replacement = "NMS.2",strk$Sample)

table(strk$Sample)

todel <- table(names(strk))[table(names(strk))>=2]
strk <- strk[!names(strk)%in%names(todel)]
strk <- strk[strk$AAChange.refGene !="UNKNOWN"]
strk <- strk[strk$Gene.refGene != "MUC3A"]
strk <- strk[strk$Gene.refGene != "MUC5AC"]
strk[strk$Gene.refGene=="MUC5AC"]
strk[strk$Gene.refGene=="MUC3A"]
strk <- strk[!grepl("PRAMEF",strk$Gene.refGene)]
strk <- strk[!grepl("LILR",strk$Gene.refGene)]
strk <- strk[!strk$Gene.refGene %in% c("OR52E5","OR52L1","SMPD1")]

message("SNVs with VAF < ", 0.10, " ", paste0(names(table(strk$Sample[strk$VAF<0.1])), ":", table(strk$Sample[strk$VAF<0.1]),sep = " "))
message("SNVs with VAF < ", 0.15, " ", paste0(names(table(strk$Sample[strk$VAF<0.15])), ":", table(strk$Sample[strk$VAF<0.15]),sep = " "))
message("SNVs with VAF < ", 0.20, " ", paste0(names(table(strk$Sample[strk$VAF<0.20])), ":", table(strk$Sample[strk$VAF<0.20]),sep = " "))


hist(strk$VAF[strk$VAF>0.10])

#We filter those variants with low VAF. By this way we are selecting the variants of the originating clone.
strk <- strk[which(strk$VAF > min.VAF)]


message("SNVs after variants with VAF < ", min.VAF, " ", paste0(names(table(strk$Sample)), ":", table(strk$Sample), sep = " "))

table(strk$Sample)
todel <- which(strk$ExonicFunc.refGene =="nonframeshift_deletion" & strk$avsnp150 != ".")
todel <- c(todel,which(strk$ExonicFunc.refGene =="nonframeshift_insertion" & strk$avsnp150 != "."))

strk <- strk[-todel]
mut <- table(strk$Sample)
barplot(mut, las =2)
table(strk$Sample)

#Extracting REF and ALT alleles.
ref.allele <- names(strk)
ref.allele <- strsplit(x= ref.allele,"_")
ref.allele <- lapply(ref.allele,function(x) x[2])
ref.allele <- lapply(ref.allele,function(x) strsplit(x,"/"))
alt.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [2])))
ref.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [1])))
strk$Ref <- ref.allele
strk$Alt <- alt.allele


strelka.variants <- split(strk,strk$Sample)
# save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.3.RData"))
load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.3.RData"))
# save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_WES.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.3.RData"))


lapply(strelka.variants,length)

# # save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.2.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.2.RData"))

colnames(mcols(strelka.variants$`90-8TL`))
mut.signature.df <- data.frame(CHROM = as.character(seqnames(strk)),
                               POS = start(strk),
                               REF = strk$REF, ALT = strk$Alt,
                               SAMLEID = strk$Sample) 
# write.table(mut.signature.df,file = file.path(execution.dir,"results","STS-26T","Strelka_GermVariants2.9.10","results",paste0("STS-26T_MutSignVariants.csv")),sep = "\t", col.names = T, row.names = F)
write.table(mut.signature.df,file = file.path(execution.dir,"results",paste0("AllSamples_MutSignVariants.3.csv")),sep = "\t", col.names = T, row.names = F)
# write.table(mut.signature.df,file = file.path(execution.dir,"results",paste0("AllSamples_MutSignVariants.TMB.csv")),sep = "\t", col.names = T, row.names = F)
# 
# 
# ## Tumor Mutation Burden
# strk.tmb <- strk
# strk.tmb <- strk.tmb[strk.tmb$Func.refGene == "exonic"]
# strk.tmb <- strk.tmb[strk.tmb$ExonicFunc.refGene %in% c("frameshift_insertion","frameshift_deletion", "nonsynonymous_SNV", "nonframeshift_deletion","nonframeshift_insertion", "stopgain","starloss")]
# # strk.tmb <- strk.tmb[strk.tmb$avsnp150 == "."]
# table(duplicated(names(strk.tmb)))
# mut <- table(strk.tmb$Sample)
# tmb.mb <- table(strk.tmb$Sample)
# 
# tmb.samp <- strk.tmb[strk.tmb$Sample =="90-8TL"]
# range <- toGRanges("chr1:1-1000000")
# subsetByOverlaps(tmb.samp,range)
# n <- 1e6
# gr <- GRanges(seqnames = as.character(seqlevels(strk.tmb)),))
# tmb.mb/3000
# barplot(tmb.mb, las =2)
# strk <- strk.tmb





# # save(strelka.files, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_NewFilters.RData"))
# # load(file.path(execution.dir, "results/All_samples_Strelka_Variants_NewFilters.RData"))
# 
# # load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.2.RData"))
# 
# #save specific variants
# for(i in seq_len(length(sample.names))){
#   sn <- sample.names[i]
#   strk <- strelka.variants[[sn]]
#   # strk <- strelka.files[[sn]]
#   rowRanges.celltype.df <- toDataframe(strk)
#   rowRanges.celltype.df$Variants <- rownames(rowRanges.celltype.df)
#   # rowRanges.celltype.df<- rowRanges.celltype.df[,c(1:26,28)]
#   write.table(rowRanges.celltype.df, file.path(execution.dir,"results",sn,"Strelka_GermVariants2.9.10/results", paste0(sn,"_ExonicsnpAlteration_toKeepFiltering.2.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   
#   
#   rowRanges.celltypes <- subsetByOverlaps(strk,gene.markers.gf)
#   rowRanges.celltype.df <- toDataframe(rowRanges.celltypes)
#   rowRanges.celltype.df$Variants <- rownames(rowRanges.celltype.df)
#   # rowRanges.celltype.df<- rowRanges.celltype.df[,c(1:26,28)]
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#   write.table(rowRanges.celltype.df, file.path(execution.dir,"results",sn,"Strelka_GermVariants2.9.10/results", paste0(sn,"_ExonicsnpAlteration_specificGenes.toKeepFiltering.2.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   
#   rowRanges.celltypes <- subsetByOverlaps(strk,all.genes.cosmic.gr)
#   rowRanges.celltype.df <- toDataframe(rowRanges.celltypes)
#   rowRanges.celltype.df$Variants <- rownames(rowRanges.celltype.df)
#   # rowRanges.celltype.df<- rowRanges.celltype.df[,c(1:26,28)]
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#   write.table(rowRanges.celltype.df, file.path(execution.dir,"results",sn,"Strelka_GermVariants2.9.10/results", paste0(sn,"_ExonicsnpAlteration_cosmicGenes.toKeepFiltering.2.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   
# }
# 
# 
# 
# # 
# # 
# # 
# # colnames(rowRanges.celltype.df)
# #   rowRanges.celltype.df <- toDataframe(rowRanges.celltypes)
# #   
# 
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_NewFiltering.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_toKeepFiltering.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   
#   ## Selecting specific genes
#   rowRanges.celltypes <- subsetByOverlaps(rowRanges.celltypes,gene.markers.gf)
#   rowRanges.celltype.df <- toDataframe(rowRanges.celltypes)
#   rowRanges.celltype.df$Variants <- rownames(rowRanges.celltype.df)
#   rowRanges.celltype.df<- rowRanges.celltype.df[,c(1:27)]
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#   write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.toKeepFiltering.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
#   
#   # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(cell.line,"_Alteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
  
  
  # 
  # 
  # rowRanges.celltypes$Variant_sel_status[patho >=5] <- "Pathogenic"
  # rowRanges.celltypes$Variant_sel_status[rowRanges.celltypes$Func.refGene =="splicing"] <- "Pathogenic"
  # rowRanges.celltypes$Variant_sel_status[rowRanges.celltypes$CLNSIG %in%c("Pathogenic","Likely_pathogenic")] <- "Pathogenic"
  # 
  # tokeep <- rep(FALSE,length(rowRanges.celltypes))
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # #We also get those variants no annotated in any predictor as "."
  # uncert <- rowSums(df.pred == ".")
  # rowRanges.celltypes$Variant_sel_status[uncert ==8] <- "."
  # 
  # sum.pat.uncert <-patho+uncert
  # rowRanges.celltypes$Variant_sel_status[sum.pat.uncert ==8] <- "."
  # 
  # rowRanges.celltypes$Variant_sel_status[sum.pat.uncert <5] <- ""
  # 
  # 
  # #we select those pathogenic and "." variants
  # rowRanges.celltypes <- rowRanges.celltypes[rowRanges.celltypes$Variant_sel_status%in% c(".","Pathogenic")]
  # 
  # seqlevels(rowRanges.celltypes) <- sortSeqlevels(seqlevels(rowRanges.celltypes))
  # rowRanges.celltypes <- sort(rowRanges.celltypes)
  # strelka.files[[sample]] <-rowRanges.celltypes
  # 
  # rowRanges.celltype.df <- toDataframe(rowRanges.celltypes)
  # rowRanges.celltype.df$Variants <- rownames(rowRanges.celltype.df)
  # rowRanges.celltype.df<- rowRanges.celltype.df[,c(27,1:26)]
  # 
  # # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  # # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_NewFiltering.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  # 
  # ## Selecting specific genes
  # rowRanges.celltypes <- subsetByOverlaps(rowRanges.celltypes,gene.markers.gf)
  # rowRanges.celltype.df <- toDataframe(rowRanges.celltypes)
  # rowRanges.celltype.df$Variants <- rownames(rowRanges.celltype.df)
  # rowRanges.celltype.df<- rowRanges.celltype.df[,c(27,1:26)]
  # # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration_specificGenes.NewFilters.csv")), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  # 
  # # write.table(rowRanges.celltype.df, file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(cell.line,"_Alteration_specificGenes.2.csv")), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  # 
  



#recurrence of SNVs
# 
# recurrence.list <- list()
# gr <- c()
# for(i in seq_len(length(sample.names))){
#   sample <- sample.names[i]
#   message("Getting Setrelka variants of specific Genes from ", sample)
#   # Directories
#   result.dir <- file.path(execution.dir,"results",sample)
#   
#   df <-read.table(file.path(result.dir,"Strelka_GermVariants2.9.10/results", paste0(sample,"_ExonicsnpAlteration.csv")), sep = "\t", header = TRUE)
#   df$Sample <- sample
#   gr <- c(gr,toGRanges(df))
# }

