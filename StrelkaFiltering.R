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

#################################### Cell Lines ############################
# Directory
execution.dir <-"/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS"
# execution.dir <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WES"

sample.data <- read.table(file = file.path(execution.dir,"Sample.info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample.names <- sample.data$Sample.name
names(sample.names) <- sample.names
prop.thr <- 0.7

########################### Tumors #################################
# Directory
# execution.dir <-"/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WGS"

# sample.names <- read.table(file = file.path(execution.dir,"sample.info.tumors.csv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# sample.names <- sample.names$Sample.name[c(2,3,8,12,14,17,22,23)]
# sample.names <- sample.names$Sample.name[c(2,3,8,12)]
# prop.thr <- 0


# Reference genome used
ref.genome.fq <- "/imppc/labs/eslab/mmagallon/Genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
ref.genome <- "hg38"

#repeatMasker 
# rptmskr <- read.table(file.path("/imppc/labs/eslab/mmagallon/Annotatiions/RepeatMasker/rmskr_hg38.bed"),header = FALSE, stringsAsFactors = FALSE)
# rptmskr <- toGRanges(rptmskr, genome = ref.genome)

gene.markers <- read.table(file.path("./special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
gene.markers.gf <- toGRanges(gene.markers$V2, genome ="hg38")
mcols(gene.markers.gf) <- gene.markers$V1
colnames(mcols(gene.markers.gf)) <- "Genes"
gene.markers.gf <- sort(gene.markers.gf)
regions.selected <- gene.markers.gf + 1e6

# TSG  COSMIC list: this was downloaded from COSMIC website but make reference to uk.sanger: http://cancer.sanger.ac.uk/census,
# this is de webpage to download this information 
# https://cog.sanger.ac.uk/cosmic/GRCh38/cosmic/v95/cancer_gene_census.csv?AWSAccessKeyId=KRV7P7QR9DL41J9EWGA2&Expires=1642007506&Signature=oVGW10M4sf3rGTeU06M1bSFFlys%3D

all.genes.cosmic <- read.table(file.path("./Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]


### Loading and filtering the strelka files

#Parameters to filte strelka output

info.param <- c("Gene.refGene",
                "ExonicFunc.refGene",
                "Func.refGene",
                "avsnp150",
                "AAChange.refGene",
                "CLNALLELEID",
                "CLNDN",
                "CLNDISDB",
                "CLNREVSTAT",
                "CLNSIG",
                "SIFT_pred",
                "Polyphen2_HDIV_pred",
                "Polyphen2_HVAR_pred",
                "LRT_pred",
                "FATHMM_pred",
                "MutationTaster_pred",
                "MutationAssessor_pred",
                "InterVar_automated",
                "AF",
                "AF_popmax","GeneDetail.refGene")
geno.param <- "AD"


#Min VAF to filter
min.VAF <- 0.1

strelka.files <- list()

for(sn in seq_len(length(sample.names))){
  sample <- sample.names[sn]
  message("Getting Setrelka variants of specific Genes from ", sample)
  # Directories
  result.dir <- file.path(execution.dir,"results",sample)
  strelka.res.path <- file.path(result.dir, "Strelka_GermVariants2.9.10/results/variants", paste0(sample, "_annotated_PASS.vcf.gz"))
  
  # strelka.res.path <- file.path(result.dir, "Strelka_GermVariants2.9.10/results/variants", paste0(sample, "_annotated.vcf.gz"))
  
  # We select only those variants present in exons or affecting splicing and we delete the synonymous_SNVs
  filt <- FilterRules(list(
    Func.refGene = function(x) as.character(unlist(info(x)$Func.refGene)) %in% c("splicing","exonic"),
    ExonicFunc.refGene = function(x) as.character(unlist(info(x)$ExonicFunc.refGene))%in% c(".","frameshift_insertion","frameshift_deletion", "stopgain", "starloss","nonsynonymous_SNV","nonframeshift_deletion","nonframeshift_insertion")))
  
  
  filt2 <- filterVcf(strelka.res.path, "hg38", tempfile() , param =ScanVcfParam(fixed = NA, geno = geno.param, info=info.param) ,filters=filt, verbose = FALSE)
  
  #Loading strelka variants
  vcf <- VariantAnnotation::readVcf(file = filt2, genome = ref.genome, param = ScanVcfParam(fixed = NA, geno =geno.param,info=info.param ))
  
  #Remove temp file filt2 
  file.remove(filt2)
  
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
  
  # Starting filtering. Selecting AF< 0.01
  message("Starting filtering. Selecting AF< 0.01")
  
  rowRanges.celltype$AF_popmax <-   as.numeric(rowRanges.celltype$AF_popmax)
  rowRanges.celltype$AF_popmax[is.na(rowRanges.celltype$AF_popmax)] <- "."
  rowRanges.celltypes <- rowRanges.celltype[rowRanges.celltype$AF_popmax<0.01 | rowRanges.celltype$AF_popmax =="."]
  
  message("SNVs after AF filtering ", length(rowRanges.celltypes))
  
  # We fileter benign or likely benign variants in CLNSIG
  ll <- lapply(rowRanges.celltypes$CLNSIG,function(x)length(x))
  rowRanges.celltypes$CLNSIG[ll >=2 ]<-lapply(rowRanges.celltypes$CLNSIG[ll>=2],function(x)paste0(x,collapse = ","))
  rowRanges.celltypes <- rowRanges.celltypes[!grepl("benign",rowRanges.celltypes$CLNSIG, ignore.case = T)]
  
  message("SNVs after CLNSIG benign filtering ", length(rowRanges.celltypes))
  
  # We delete those benign or likely_benign variants annotated in InterVar_automated
  rowRanges.celltypes <- rowRanges.celltypes[!grepl("benign",rowRanges.celltypes$InterVar_automated,ignore.case = T)]
  
  message("SNVs after InterVar_automated benign or likely benign filtering ", length(rowRanges.celltypes))
  
  # We delete those variants with UNKNOWN in AAChange.refGene
  rowRanges.celltypes <- rowRanges.celltypes[rowRanges.celltypes$AAChange.refGene !="UNKNOWN"]
  message("SNVs after filtering UNKNOWN variants ", length(rowRanges.celltypes))
  
  
  #Then, we decided a variant will be pathogenic if more than 5 predictors say that.
  df <- data.frame(ExonicFunc.refGene = as.character(rowRanges.celltypes$ExonicFunc.refGene), 
                   Func.refGene = as.character(rowRanges.celltypes$Func.refGene) , 
                   SIFT_pred = as.character(rowRanges.celltypes$SIFT_pred),
                   Polyphen2_HDIV_pred = as.character(rowRanges.celltypes$Polyphen2_HDIV_pred),
                   Polyphen2_HVAR_pred = as.character(rowRanges.celltypes$Polyphen2_HVAR_pred),
                   LRT_pred = as.character(rowRanges.celltypes$LRT_pred),
                   FATHMM_pred = as.character(rowRanges.celltypes$FATHMM_pred),
                   MutationTaster_pred = as.character(rowRanges.celltypes$MutationTaster_pred),
                   MutationAssessor_pred = as.character(rowRanges.celltypes$MutationAssessor_pred))
  # InterVar_automated = as.character(rowRanges.celltypes$InterVar_automated))
  
  df$SIFT_pred[df$SIFT_pred== "D"] <- "Pathogenic"
  df$Polyphen2_HDIV_pred[df$Polyphen2_HDIV_pred%in%c("D","P")] <- "Pathogenic"
  df$Polyphen2_HVAR_pred[df$Polyphen2_HVAR_pred%in%c("D","P")] <- "Pathogenic"
  df$LRT_pred[df$LRT_pred == "D"] <- "Pathogenic"
  df$MutationTaster_pred[df$MutationTaster_pred %in% c("A","D","P")] <- "Pathogenic"
  df$MutationAssessor_pred[df$MutationAssessor_pred %in% c("H","M")] <- "Pathogenic"
  df$FATHMM_pred[df$FATHMM_pred == "D"] <- "Pathogenic"
  # df$InterVar_automated[df$InterVar_automated == "Likely_pathogenic"] <- "Pathogenic"
  
  #Counting the number of pathogenic predictors
  df.pred <- df[,-c(1,2)]
  patho <- rowSums(df.pred == "Pathogenic") 
  uncert <- rowSums(df.pred == ".")
  prop <- (patho)/(7-uncert)
  # prop[584]
  #NaN results are normaly population variants but we filter common variants later.
  #now we keep those NaN values.
  rowRanges.celltypes$prop_pathogenic_pred <- prop
  tokeep <- rep(FALSE,length(rowRanges.celltypes))
  
  tokeep[rowRanges.celltypes$Func.refGene == "splicing"] <- TRUE
  tokeep[prop>prop.thr] <- TRUE
  
  tokeep[is.nan(prop)] <- TRUE
  tokeep[grepl(pattern = "pathogenic", rowRanges.celltypes$InterVar_automated, ignore.case = T)] <- TRUE
  tokeep[grepl(pattern = "pathogenic", rowRanges.celltypes$CLNSIG,ignore.case = T)] <- TRUE
  rowRanges.celltypes <- rowRanges.celltypes[tokeep]
  
  message("SNVs after filtering  variants < ", prop.thr," for pathogenic predictors ", length(rowRanges.celltypes))
  
  seqlevels(rowRanges.celltypes) <- sortSeqlevels(seqlevels(rowRanges.celltypes))
  rowRanges.celltypes <- sort(rowRanges.celltypes)
  rowRanges.celltypes <- filterChromosomes(rowRanges.celltypes)
  
  strelka.files[[sample]] <-rowRanges.celltypes
}
# save(strelka.files, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFileter.3.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFileter.2.RData"))

# After filtering the variants per sample, now we proceed to filter common variants
#If data is comming from cell lines execute the following three lines and follow the other lines down
strk <- unlist(as(strelka.files,"GRangesList"))
names(strk)<-gsub(pattern = "NMS.2",replacement = "NMS-2",names(strk))
names(strk)<-gsub(pattern = "SNF96.2",replacement = "SNF96-2",names(strk))

# Else omit the three before lines and start from here
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


hist(strk$VAF[strk$VAF>0.20])

#We filter those variants with low VAF. By this way we are selecting the variants of the originating clone.
strk <- strk[which(strk$VAF > min.VAF)]


todel <- which(strk$ExonicFunc.refGene =="nonframeshift_deletion" | strk$ExonicFunc.refGene =="nonframeshift_insertion")
if(length(unique(strk$avsnp150[todel])) !=1){
  todel <- todel[strk$avsnp150[todel] != "."]
  strk <- strk[-todel]
  
}else{
  strk <- strk
}

mut <- table(strk$Sample)
barplot(mut, las =2)
ta(strk)
tt <- strk[strk$Sample=="STS-26T"]
tt[tt$Gene.refGene=="CDKN2A"]

# Extracting REF and ALT alleles.
ref.allele <- names(strk)
ref.allele <- strsplit(x= ref.allele,"_")
ref.allele <- lapply(ref.allele,function(x) x[2])
ref.allele <- lapply(ref.allele,function(x) strsplit(x,"/"))
alt.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [2])))
ref.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [1])))
strk$Ref <- ref.allele
strk$Alt <- alt.allele


strelka.variants <- split(strk,strk$Sample)
length(strelka.variants$`STS-26T`)
strelka.variants


save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.3.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.RData"))

load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.3.RData"))
lapply(strelka.variants,length)

# save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_WES_toKeepFilters.2.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_WES_toKeepFilters.2.RData"))

######## Selecting those common variants present in both WGS and WES data  #################
strk.lines.wgs <- strelka.variants

# saving specific genes mutations and cosmic genes mutations
for(i in seq_len(length(strk.lines.wgs))){
  sn <- names(strk.lines.wgs)[i]
  wgs <- strk.lines.wgs[[sn]]
  wgs.df <- toDataframe(wgs)
  wgs.df$Variants <- rownames(wgs.df)
  write.table(wgs.df, file.path(execution.dir,"results", sn, "Strelka_GermVariants2.9.10/results/",paste0( sn, "_ExonicsnpAlteration.toKeepFiltering.3.csv")), row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  specific.genes <- toDataframe(subsetByOverlaps(wgs,gene.markers.gf))
  specific.genes$Variants <- rownames(specific.genes)
  write.table(specific.genes, file.path(execution.dir,"results", sn, "Strelka_GermVariants2.9.10/results/",paste0( sn, "_ExonicsnpAlteration_specificGenes.toKeepFiltering.3.csv")), row.names = FALSE, col.names = TRUE, quote = FALSE)
  cosmic.genes <- toDataframe(subsetByOverlaps(wgs, all.genes.cosmic.gr))
  cosmic.genes$Variants <- rownames(cosmic.genes)
  
  write.table(cosmic.genes, file.path(execution.dir,"results", sn, "Strelka_GermVariants2.9.10/results/",paste0( sn, "_ExonicsnpAlteration_cosmicGenes.toKeepFiltering.3.csv")), row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}

