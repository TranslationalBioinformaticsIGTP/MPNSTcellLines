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
# source(file = "/imppc/labs/eslab/mmagallon/Projects/cliffhunter/SV.functions.R")
source(file = "./Integrative_Biology_Functions.R")

#################### Parameters & Directories #####################

############## Cell Lines ######
#Directory
execution.dir <-"./MPNST_cellLines/WGS"

sample.data <- read.table(file = file.path(execution.dir,"Sample.info.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sample.names <- sample.data$Sample.name
names(sample.names) <- sample.names
num.rep.samples <- 1
wes.wgs <- "wgs"

#Genome
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
                "AF",
                "AF_popmax","GeneDetail.refGene")


geno.param <- c("AD","DP")


#Min VAF to filter
min.VAF <- 0.1

# loading somatic rs present in icgc and cosmic
all.rs <- read.table(file.path("/imppc/labs/eslab/mmagallon/Annotatiions/all.rsnames.icgcAndCosmic.txt"), header = F)

strelka.files <- list()

for(sn in seq_len(length(sample.names))){
  sample <- sample.names[sn]
  message("Getting Setrelka variants of specific Genes from ", sample)
  # Directories
  result.dir <- file.path(execution.dir,"results",sample)
  strelka.res.path <- file.path(result.dir, "Strelka_GermVariants2.9.10/results/variants", paste0(sample, "_annotated_PASS.vcf.gz.tbi"))
  
  rowRanges.celltypes <- loadingAnnotStrelkaVCFforMutSig(vcf.annotated.file = strelka.res.path, 
                                            genome = ref.genome,
                                            info.param = info.param,
                                            geno.param = geno.param,
                                            wes.wgs = wes.wgs,
                                            add.vaf = TRUE,
                                            verbose = TRUE) 

   #Filtering per sample
  message("Starting filtering. Selecting AF< 0.01")
  
  rowRanges.celltypes <- rowRanges.celltypes[rowRanges.celltypes$AF_popmax<0.01 | rowRanges.celltypes$AF_popmax =="."]
  
  message("SNVs after AF filtering ", length(rowRanges.celltypes))
  
  # We delete those variants with UNKNOWN in AAChange.refGene
  rowRanges.celltypes <- rowRanges.celltypes[rowRanges.celltypes$AAChange.refGene !="UNKNOWN"]
  message("SNVs after UNKNOWN filtering ", length(rowRanges.celltypes))
  
  
  seqlevels(rowRanges.celltypes) <- sortSeqlevels(seqlevels(rowRanges.celltypes))
  rowRanges.celltypes <- sort(rowRanges.celltypes)
  rowRanges.celltypes <- filterChromosomes(rowRanges.celltypes)
  rowRanges.celltypes$Sample <- sample
  rowRanges.celltypes$Variant <- names(rowRanges.celltypes)
  
  strelka.files[[sample]] <-rowRanges.celltypes
}


# save(strelka.files, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignature.3.RData"))#Paper
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignature.3.RData"))

save(strelka.files, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignature.v4.RData"))#function
load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignature.v4.RData"))

# ss <- data.frame(strelka.files$S462)

# After filtering the variants per sample, now we proceed to filter common variants
strk <- unlist(as(strelka.files,"GRangesList"))
names(strk) <-strk$Variant

table(strk$Sample)
message("Total SNPs ", length(strk))

# todel <- table(names(strk))[table(names(strk))>=num.rep.samples]
todel <- names(table(names(strk)))[table(names(strk))>num.rep.samples]
strk <- strk[!names(strk) %in% todel]
table(strk$Sample)
message(length(strk)," SNV initials")



#Deletion of variants present in highly variable genes
strk <- strk[strk$Gene.refGene != "MUC3A"]
strk <- strk[strk$Gene.refGene != "MUC5AC"]
strk <- strk[!grepl("PRAMEF",strk$Gene.refGene)]
strk <- strk[!grepl("LILR",strk$Gene.refGene)]
strk <- strk[!strk$Gene.refGene %in% c("OR52E5","OR52L1","SMPD1")]
table(strk$Sample)
message("SNVs after UNKNOWN filtering ", length(strk))

#We filter those variants with low VAF. By this way we are selecting the variants of the originating clone.
message("SNVs with VAF < ", 0.10, " ", paste0(names(table(strk$Sample[strk$VAF<0.1])), ":", table(strk$Sample[strk$VAF<0.1]),sep = " "))
message("SNVs with VAF < ", 0.15, " ", paste0(names(table(strk$Sample[strk$VAF<0.15])), ":", table(strk$Sample[strk$VAF<0.15]),sep = " "))
message("SNVs with VAF < ", 0.20, " ", paste0(names(table(strk$Sample[strk$VAF<0.20])), ":", table(strk$Sample[strk$VAF<0.20]),sep = " "))

tokeep <- rep(FALSE,length(strk))
tokeep[which(strk$VAF > min.VAF)] <- TRUE
tokeep[is.na(strk$VAF)] <- TRUE
strk <- strk[tokeep]
table(strk$Sample)
message("SNVs after variants with VAF < ", min.VAF, " ", paste0(names(table(strk$Sample)), ":", table(strk$Sample), sep = " "))

todel <- which(strk$ExonicFunc.refGene =="nonframeshift_deletion" & strk$avsnp150 != ".")
todel <- c(todel,which(strk$ExonicFunc.refGene =="nonframeshift_insertion" & strk$avsnp150 != "."))
strk <- strk[-todel]

mut <- table(strk$Sample)
barplot(mut, las =2)
table(strk$Sample)

#Filtering out those variant with rs not in icgc nor cosmic.
known <- which(strk$avsnp150!=".")
som.rs <- strk$avsnp150[known][strk$avsnp150[known]%in% all.rs$V1]
strk <- strk[strk$avsnp150 %in% c( som.rs, ".")]

mut <- table(strk$Sample)
barplot(mut,las=2)

strelka.variants <- split(strk, strk$Sample)
# save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.3.RData"))#paper
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.3.RData"))
save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.v4.RData"))#Filter icgc
load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.v4.RData"))

# save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_WES.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_MutSignatureFiltered.3.RData"))


unlist(lapply(strelka.variants,length))

# # save(strelka.variants, file = file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.2.RData"))
# load(file.path(execution.dir, "results/All_samples_Strelka_Variants_toKeepFilters.2.RData"))

mut.signatures.df <- data.frame(CHROM = as.character(seqnames(strk)),
                               POS.START = start(strk),
                               POS.END = end(strk),
                               REF = strk$REF, 
                               ALT = strk$ALT,
                               SAMPLEID = strk$Sample) 
# write.table(mut.signature.df,file = file.path(execution.dir,"results","STS-26T","Strelka_GermVariants2.9.10","results",paste0("STS-26T_MutSignVariants.csv")),sep = "\t", col.names = T, row.names = F)
# write.table(mut.signature.df,file = file.path(execution.dir,"results",paste0("AllSamples_MutSignVariants.3.csv")),sep = "\t", col.names = T, row.names = F)#paper
write.table(mut.signatures.df,file = file.path(execution.dir,"results",paste0("AllSamples_MutSignVariants.v4.csv")),sep = "\t", col.names = T, row.names = F)#functions

