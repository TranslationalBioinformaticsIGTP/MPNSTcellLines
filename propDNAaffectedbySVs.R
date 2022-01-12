# Packages Needed
library(CopyNumberPlots)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(regioneR)

######################## Functions #########################################
source(file = "/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/loadingLumpySVs.R")
source(file = "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS/funtionsPaperCellLines.R")
source(file= "/imppc/labs/eslab/mmagallon/Projects/RNA-Seq-timecourse.2/Analysis/rna_seq_Functions.R")
################## Parameters #######################
execution.dir <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS"

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
#loading data
load(file= file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.v3.RData"))
load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.v4.RData")))
filt.list.icr.inv <- lapply(filt.list, function(x) x[["AllRegionsFiltered"]])
filt.list.icr.inv <- lapply(filt.list.icr.inv, function(x) x[x$SVTYPE %in% c("INV","ICR")])

for(i in seq_len(length(sample.names))){
  sn <- sample.names[i]
  cliff <- join.cliff.lumpy[[sn]]
  filt.list.icr.inv[[sn]] <- unique(c(filt.list.icr.inv[[sn]],cliff))
}
filt.list.icr.inv <- lapply(filt.list.icr.inv, function(x) x[x$SVTYPE %in% c("INV","ICR")])

#cantidad de genoma afectada por SVs
extend1 <- 1e6
svs.mb <- list()
for(i in seq_len(length(sample.names))){
  sn <- sample.names[i]
  ss <- extendRegions(trim(filt.list.icr.inv[[sn]]),extend.start =extend1, extend.end = extend1)
  svs.mb[[sn]] <- reduce(ss)
}

res <- unlist(lapply(svs.mb,function(x){sum(width(x))}))/3e9
round(res,digits = 3)
barplot(res, las =2)
dev.off()
sn <- "SNF96.2"
for(i in seq_len(length(sample.names))){
  sn <- sample.names[i]
  sv.sn <- filterChromosomes(svs.mb[[sn]])
  sv.sn <- sv.sn
  sv.sn.subset <- subsetByOverlaps(filt.list.icr.inv[[sn]], sv.sn)
  chr.end <- sv.sn.subset$END
  chrs <- unique(c(chr,unlist(lapply(strsplit(chr.end,":"),function(x) x[1]))))
  chrs <- seqlevels(sv.sn)[seqlevels(sv.sn)%in%chrs]
  kp <-plotKaryotype(chromosome = chrs, genome = genome, main = sn) 
  kpPlotLinks(kp, data = sv.sn.subset, data2 = toGRanges(sv.sn.subset$END))
  kpPlotRegions(kp, sv.sn, col = transparent("grey77"))
  
}

sv.sn.subset[seqnames(sv.sn.subset) == "chr9"]
sv.sn <- svs.mb[[sn]]
sv.sn.subset <- subsetByOverlaps(filt.list.icr.inv[[sn]], sv.sn)
chr.end <- sv.sn.subset$END
chrs <- unique(c(chr,unlist(lapply(strsplit(chr.end,":"),function(x) x[1]))))

length(sv.sn.subset)
kp <-plotKaryotype(chromosome = chrs, genome = genome) 
kpPlotRegions(kp, sv.sn, col = transparent("grey77"))
kpPlotLinks(kp, data = sv.sn.subset, data2 = toGRanges(sv.sn.subset$END))

