##################################################
#     Script Figures MPNST cell lines paper      #
##################################################

############### Packages needed  #################
# Packages Needed
library(CopyNumberPlots)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38.masked)

######################## Functions #########################################
source(file = "./loadingLumpySVs.R")
source(file = "./funtionsPaperCellLines.R")
################## Parameters #######################
execution.dir <- "./MPNST_cellLines/WGS"

#Cosmic and Gene files
#Selection of interesting genes
gene.markers <- read.table(file.path("./special.genes.txt"), header = FALSE, sep = " ", stringsAsFactors = FALSE)
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
all.genes.cosmic <- read.table(file.path("./Census_allThu Sep 16 11_23_32 2021.tsv"), sep = "\t",header = T)
all.genes.cosmic$Genome.Location <-  paste0("chr",all.genes.cosmic$Genome.Location)
all.genes.cosmic <- all.genes.cosmic[nchar(all.genes.cosmic$Genome.Location)>7,]
all.genes.cosmic.gr <- toGRanges(all.genes.cosmic$Genome.Location)
mcols(all.genes.cosmic.gr) <- all.genes.cosmic[,c(1,2,5,6:20)]

all.genes.cosmic[545,]
tumor.sup.gr <- all.genes.cosmic.gr[grepl("TSG", all.genes.cosmic.gr$Role.in.Cancer)]



#Fusion Genes Cosmic

table(all.fusion.genes.cosmic$PRIMARY_HISTOLOGY)
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


# gene.markers <- toGRanges(data.frame(chr= c("chr9", "chr17", "chr17", "chr17", "chr11"), start = c(21967753,31094927, 31937007,7668230, 86245050), end = c(21975098, 31377677, 32001038, 7687366,86278810), symbol = c("CDKN2A", "NF1", "SUZ12", "TP53", "EED")))
# gene.markers <- toGRanges(gene.markers.df$V2, genome="hg38")
# mcols(gene.markers) <- gene.markers.df$V1
# colnames(mcols(gene.markers)) <- "symbol"
# gene.markers <- sort(gene.markers)





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

######################### Loading the data Needed ##############################

###CNVkit files
#CNVkit lrr
load(file= file.path(execution.dir, "results/CNVkit_lrr_AllSamples.RData"))
#CNVkit CN
load(file= file.path(execution.dir, "results/CNVkit_cn_AllSamples.RData"))

### SNParray 

##### loading GAP SNP-array data from MPNST Cell lines
data.dir <- ""
result.dir <- "./MPNST_cellLines/SNP-array/results/paperCDKN2A"

cell.lines <- read.table("./MPNST_cellLines/SNP-array/cellLinesDataBernatDir.csv",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
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


# We have to load GAP snp-array data from Tumors
#Loading segment data
sample.info.tumors <- read.table(file.path("./SNP_ARRAY_MPNST.csv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
file.path.mpnst <- sample.info.tumors$File_Path
data.path.mpnst <- file.path(sample.info.tumors$Data_path,"data")

off.name <- sample.info.tumors$Official_Name
sample.names.tumor <- sample.info.tumors$Sample_Name
names(file.path.mpnst) <- sample.names.tumor
names(off.name)<- sample.names.tumor
names(data.path.mpnst) <- sample.names.tumor

data.dir <- data.path.mpnst
result.dir <- "./SNPArrays/results/hg38SnpArray/"
sample.names

segment.files.tumor <- list()
snp.files.tumor <- list()
i=4
for (i in seq_len(length(sample.names.tumor))){
  #segment data 
  sn <- sample.names.tumor[i]
  segment.file <- file.path.mpnst[sn]
  off.nm <- off.name[sn]
  
  
  if(file.exists(segment.file)){
    segment.files.tumor[[off.nm]] <- loadCopyNumberCalls(cnv.data = segment.file) 
    segment.files.tumor[[off.nm]] <- UCSCStyle(segment.files.tumor[[off.nm]])
    segment.files.tumor[[off.nm]]$cn[segment.files.tumor[[off.nm]]$cn %in% 100] <- 0
    
  } else{
    stop(segment.file, "doesn't exist")
  }
  
  #snpData
  data.dir <- data.path.mpnst[sn]
  data.dir <- file.path(data.dir,list.files(data.dir))
  snp.files.tumor[[off.nm]] <- loadSNPData(snp.data = data.dir) 
  snp.files.tumor[[off.nm]] <- UCSCStyle(snp.files.tumor[[off.nm]])
}
i=1
# Fixing GAP error in LOH
for(i in seq_len(length(snp.files.tumor))){
  segment.files.tumor[[i]]$loh[segment.files.tumor[[i]]$cn ==0] <-0
}


# Loading the SNParray data from WGS
load(file.path(execution.dir,"results","MPNST_WGS_cellLines_BAF_Filtered.RData"))

### Joined CliffHunter and Lumpy data of specific genes
load(file= file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.v4.RData"))

### Loading the table to plot with the mutational information of important genes
load(file = file.path(execution.dir,"results","table.to.plot.toKeep.3.RData"))

################################ Images Paper ########################################

## F1 B
#Plot Params
cn.colors <- getCopyNumberColors()
cn.colors[3] <- "white"

pp <- getDefaultPlotParams(plot.type = 3)
plotDefaultPlotParams(plot.type = 3,plot.params = pp)
pp$topmargin <- 50
pp$bottommargin <- 20
pp$data1inmargin <- 5
pp$data2inmargin <- 5
pp$ideogramheight <- 8
pp$ideogramlateralmargin <- 0.003
pp$leftmargin <- 0.15
pp$rightmargin <- 0.1
pp$data1outmargin <- 20

#label sizes
cex.chr <- 5
cex.axis <- 4
cex.titel.axis <- 5

group <- c("SP","NF1","All")

for(j in seq_len(length(group))){
  gr <- group[j]
  if(gr =="SP"){
    segment.files.tumor.sp <- segment.files.tumor[grepl("SP", names(segment.files.tumor))]
    segment.files.sp <- segment.files[names(segment.files)%in% c("STS-26T","HS-Sch-2","HS-PSS") ]
    ymax = 4
    label.axis <- c("4","0","4")
    
  }else if (gr == "NF1"){
    segment.files.tumor.sp <- segment.files.tumor[!grepl("SP", names(segment.files.tumor))]
    segment.files.sp <- segment.files[!names(segment.files)%in% c("STS-26T","HS-Sch-2","HS-PSS") ]
    ymax = 4
    label.axis <- c("4","0","4")
    
  }else{
    segment.files.sp <- segment.files
    segment.files.tumor.sp <- segment.files.tumor
    ymax = 4
    label.axis <- c("4","0","4")
  }
  #cn of tumors and lines
  gains.tumor <- lapply(segment.files.tumor.sp, function(x){
    x[x$cn>2]
  })
  
  nt <- length(gains.tumor)
  gains.tumor <- unlist(GRangesList(gains.tumor))
  # gains.tumor$cn <- gains.tumor$cn/nt
  gain.tumors.cn <- c()
  
  for(i in seq_len(length(gains.tumor))){
    gain.tumors.cn <- c(gain.tumors.cn,rep(gains.tumor[i],(gains.tumor$cn[i]-2)))
  }
  
  gain.tumors.cn<- unlist(GRangesList(gain.tumors.cn))
  sn <- strsplit(x = names(gain.tumors.cn),"\\.")
  sn <- unlist(lapply(sn,function(x){x[1]}))
  gain.tumors.cn$sample <- sn
  
  losses.tumor <- lapply(segment.files.tumor.sp, function(x){
    x[x$cn<2]
  })
  losses.tumor <- unlist(GRangesList(losses.tumor))
  losses.tumor.cn <- c()
  i=1
  for(i in seq_len(length(losses.tumor))){
    losses.tumor.cn <- c(losses.tumor.cn,rep(losses.tumor[i],(2-losses.tumor$cn[i]))) 
  }
  losses.tumor.cn<- unlist(GRangesList(losses.tumor.cn))
  sn <- strsplit(x = names(losses.tumor.cn),"\\.")
  sn <- unlist(lapply(sn,function(x){x[1]}))
  losses.tumor.cn$sample <- sn
  
  
  #cn cell lines
  gains.cell <- lapply(segment.files.sp, function(x){
    x[x$cn>2]
  })
  nc <- length(gains.cell)
  
  gains.cell<- unlist(GRangesList(gains.cell))
  gain.cell.cn <- c()
  
  for(i in seq_len(length(gains.cell))){
    gain.cell.cn <- c(gain.cell.cn,rep(gains.cell[i],(gains.cell$cn[i]-2))) 
  }
  gain.cell.cn <- unlist(GRangesList(gain.cell.cn))
  sn <- strsplit(x = names(gain.cell.cn),"\\.")
  sn <- unlist(lapply(sn,function(x){x[1]}))
  gain.cell.cn$sample <- sn
  
  
  
  losses.cell <- lapply(segment.files.sp, function(x){
    x[x$cn<2]
  })
  
  losses.cell<- unlist(GRangesList(losses.cell))
  losses.cell.cn <- c()
  
  
  for(i in seq_len(length(losses.cell))){
    losses.cell.cn <- c(losses.cell.cn,rep(losses.cell[i],(2-losses.cell$cn[i]))) 
  }
  losses.cell.cn<- unlist(GRangesList(losses.cell.cn))
  sn <- strsplit(x = names(losses.cell.cn),"\\.")
  sn <- unlist(lapply(sn,function(x){x[1]}))
  losses.cell.cn$sample <- sn
  
  # Plot Figure 1 cell lines paper
  image.dir<-"./MPNST_cellLines/WGS/results/images.paper"
  
  
  svg(filename = file.path(image.dir, paste0("F1_celllinesAndTumors_", gr, "_alleles.2.svg")), width = 45, height = 15)
  
  kp <- plotKaryotype(genome = "hg38",chromosomes = "canonical",  plot.type = 3, plot.params = pp, cex = cex.chr, srt=45)
  # kpAddCytobandsAsLine(kp)
  kpPlotMeanCoveragebySample(kp, gain.tumors.cn, num.samples = nt, ymax = ymax,col = cn.colors[5], r0=0.5, r1=0)
  kpPlotMeanCoveragebySample(kp, gain.cell.cn, num.samples = nc, ymax = ymax, col = cn.colors[4], r0=0.5, r1=1)
  kpPlotMeanCoveragebySample(kp, losses.tumor.cn, num.samples = nt, ymax = 2, col = cn.colors[1], r0=0.5, r1=1,data.panel = 2)
  kpPlotMeanCoveragebySample(kp,losses.cell.cn, num.samples = nc, ymax = 2, col = cn.colors[2], r0=0.5, r1=0, data.panel = 2)
  
  kpAxis(kp, ymin = 0, ymax =ymax, data.panel = 1, labels = label.axis, cex = cex.axis)
  kpAbline(kp, h =c(0,0.75,1.5,2.25,3), ymax =3 )
  kpAxis(kp, ymin = 0, ymax =2, data.panel = 2, labels = c("2","0","2"), cex = cex.axis)
  kpAbline(kp, h =c(0,0.5,1,1.5,2), ymax = 2, data.panel = 2 )
  
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Gains", side = "left", cex=cex.titel.axis, r0=0, r1=1, data.panel = 1,label.margin = 0.05)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Losses", side = "left", cex=cex.titel.axis, r0=0, r1=1, data.panel = 2,label.margin = 0.05)
  kpAbline(kp, h =c(1.5),ymax = 3, lwd = 2, data.panel = 1)
  kpAbline(kp, h =c(1),ymax = 2, lwd = 2, data.panel = 2)
  
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Cell lines", side = "rigth", cex=cex.titel.axis, r0=0, r1=0.5, data.panel = 2)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Tumors", side = "rigth", cex=cex.titel.axis, r0=0.5, r1=1, data.panel = 2)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Cell lines", side = "rigth", cex=cex.titel.axis, r0=0.5, r1=1, data.panel = 1)
  karyoploteR::kpAddLabels(karyoplot = kp,labels = "Tumors", side = "rigth", cex=cex.titel.axis, r0=0, r1=0.5, data.panel = 1)
  
  karyoploteR::kpAddLabels(karyoplot = kp, labels = "# Gained alleles",side = "left", cex = cex.axis, r0=0, r1=1, data.panel = 1, srt = 90,label.margin = 0.04, pos = 1)
  karyoploteR::kpAddLabels(karyoplot = kp, labels = "# Lost alleles",side = "left", cex = cex.axis, r0=0, r1=1, data.panel = 2, srt = 90,label.margin = 0.04, pos = 1)
  
  dev.off()
  
}

####### Data Sheet and supfile figures ######

#Plot Parameters
cn.colors <- getCopyNumberColors()
cn.colors[3] <- "white"
cn.colors <- c(cn.colors, cn.colors[length(cn.colors)],cn.colors[length(cn.colors)])
names(cn.colors)[length(cn.colors)] <- 10
names(cn.colors)[length(cn.colors)-1] <- 9
cn.colors <- c(cn.colors,rep(cn.colors[length(cn.colors)],10))
names(cn.colors)[names(cn.colors)%in% "10"] <- as.character(c(10:20))
pp <- getDefaultPlotParams(plot.type = 1)
pp$topmargin <- 30
pp$bottommargin <- 20
pp$data1inmargin <- 5
pp$ideogramheight <- 2
pp$ideogramlateralmargin <- 0.008
pp$leftmargin <- 0.1

chrs <- c("canonical",unique(as.character(seqnames(sample.cn[["S462"]]))))

for(i in seq_len(length(sample.names))){
  sn <- names(sample.names)[i]
  sns <- sample.names[i]
  image.dir <- file.path(execution.dir,"results",sample.names[i],"images")

  for(chr in chrs){
    png(file.path(image.dir, paste0(sn,"_",chr,"_BAF_LRR_CN_GAP_CNVKIT_strelka_NoSexCnvkit.2.png")), width = 2500, height = 800)
    
    if(chr == "canonical"){
      kp <- plotKaryotype(plot.type = 4,   plot.params = pp, chromosomes = chr, cex= 2, genome = "hg38")
      # kpAddCytobandsAsLine(kp)
      point.cx <- 0.2
      # kpAddBaseNumbers(kp,cex=1)
    }
    
    if(chr != "canonical"){
      kp <- plotKaryotype(plot.type = 1,  main = sn, plot.params = pp, chromosomes = chr, cex= 2, genome = "hg38")
      kpAddBaseNumbers(kp, 20e6, cex = 2)
      point.cx <- 0.5
      
    }
    
    # Markers track
    
    mm <-kpPlotMarkers(kp, data = table.to.plot[[sn]], labels =table.to.plot[[sn]]$Genes ,label.color = table.to.plot[[sn]]$color,  cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1000, r0= 0, r1 = 1.4,ignore.chromosome.ends = TRUE,
                       marker.parts = c(0.9, 0.05, 0.05))
    
    cv <- mm$latest.plot$computed.values$label.position
    info <- table.to.plot[[sn]]
    info$x.pos <- NA
    for(i in seq_len(length(cv))){
      chroms <- names(cv)[i]
      cv.chr <- cv[names(cv) %in% chroms]
      cv.chr <- sort(cv.chr)
      info[seqnames(info)%in%chroms]$x.pos <-cv.chr
      
    }
    
    snv <- info[info$variant%in%c("SNV")]
    if(length(snv)>0){
      mm$beginKpPlot()
      points(x = snv$x.pos, y = rep(233,length(snv)), pch =16, col =snv$variant.col, cex = 2.5)
      mm$endKpPlot()
    }
    
    
    sv <- info[info$variant%in%c("SV")]
    if(length(sv)>0){
      mm$beginKpPlot()
      points(x = sv$x.pos, y = rep(233,length(sv)), pch =17, col =sv$variant.col, cex = 2.5)
      mm$endKpPlot()
    }
    
    #CN tracks
    at <- autotrack(current.track = 1, total.tracks = 5,margin = 0.20)
    ats <- autotrack(current.track =1, total.tracks = 2,r0 = at$r0,r1=at$r1, margin = 0.20)
    
    plotCopyNumberCalls(kp, segment.files[[sn]], r0 = ats$r0, r1=ats$r1, label.cex = 2.3, label.margin = 0.03, 
                        labels = "CN-SNParray",  loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
    
    
    ats <- autotrack(current.track = 2, total.tracks = 2,r0 = at$r0,r1=at$r1, margin = 0.20)
    plotCopyNumberCalls(kp, sample.cn[[sn]], r0 = ats$r0, r1=ats$r1, label.cex = 2.3, 
                        labels = "CN-WGS", label.margin = 0.03,  loh.height = 0, loh.column = "loh", cn.colors = cn.colors)

    #LRR
    at <- autotrack(current.track = 2:3, total.tracks = 5, margin = 0.20)
    ats <- autotrack(current.track = 1, total.tracks = 2, margin = 0.20, r0=at$r0, r1=at$r1)
    
    plotLRR(kp, snps =snp.files[[sn]], labels = "LRR-SNParray",label.cex = 2.3, label.srt = 0, label.margin = 0.06, axis.cex =1.8 , points.cex = point.cx, r0 = ats$r0, r1=ats$r1)
    
    
    ats <- autotrack(current.track = 2, total.tracks = 2, margin = 0.20, r0=at$r0, r1=at$r1)
    plotLRR(kp, snps =sample.log[[sn]], labels = "LRR-WGS",label.cex = 2.3, label.srt = 0, label.margin = 0.06, axis.cex =1.8 , points.cex = point.cx, r0 = ats$r0, r1=ats$r1,line.at.0 = FALSE)
    plotCopyNumberCallsAsLines(kp, sample.cn[[sn]], cn.column = "segment.value", style = "segments", col = "red", r0 = ats$r0, r1=ats$r1, ymin = -4, ymax=2, add.axis = FALSE, labels = NA)
    
    
    #BAF representation
    at <- autotrack(current.track = 4:5, total.tracks = 5, margin = 0.20)
    ats <- autotrack(current.track = 1, total.tracks = 2, r0=at$r0, r1=at$r1, margin = 0.20)
    
    plotBAF(kp,snps =snp.files[[sn]], baf.column = "baf",labels = "BAF-SNParray", label.cex = 2.3, label.srt = 0, label.margin = 0.06, axis.cex =1.8, points.cex = point.cx, r0 = ats$r0, r1=ats$r1,track.margin = 0.15 )
    
    ats <- autotrack(current.track = 2, total.tracks = 2, r0=at$r0, r1=at$r1)
    plotBAF(kp, snps = mpnst.snp.lists[[sns]], baf.column = "baf", labels = "BAF-WGS", label.cex = 2.3, label.srt = 0, label.margin = 0.06, axis.cex =1.8, points.cex = point.cx, r0 = ats$r0, r1=ats$r1,track.margin = 0.15 )
    
       
    dev.off()
  }
  
}

#### F2
#All celllines together
chrs="canonical"
i = 1

# Plot Parameters
cn.colors <- getCopyNumberColors()
cn.colors[3] <- "white"
cn.colors <- c(cn.colors, cn.colors[length(cn.colors)],cn.colors[length(cn.colors)])
names(cn.colors)[length(cn.colors)] <- 10
names(cn.colors)[length(cn.colors)-1] <- 9
# plotDefaultPlotParams(plot.type = 3, plot.params = pp)

pp <- getDefaultPlotParams(plot.type = 4)
pp$topmargin <- 30
pp$bottommargin <- 20
pp$data1inmargin <- 5
pp$ideogramheight <- 1
pp$ideogramlateralmargin <- 0.008
pp$leftmargin <- 0.15
pp$rightmargin <- 0.1

paper.dir <- "./MPNST_cellLines/WGS/results/images.paper"

if(!file.exists(paper.dir))dir.create(paper.dir)
sample.names <- sample.names[c(length(sample.names):1)]

chrs <- "canonical"

for(chr in chrs){
  png(file.path(paper.dir, paste0("F2_AllSamples_",chr,"_BAF_LRR_CN_GAP_CNVKIT_strelka_NoSexCnvkit.2.png")), width = 1800, height = 2000)
  
  if(chr == "canonical"){
    kp <- plotKaryotype(plot.type = 4, plot.params = pp, chromosomes = chr, cex= 1, genome = "hg38", labels.plotter = NULL)
    # kpAddCytobandsAsLine(kp)
    point.cx <- 0.1
    kpAddChromosomeNames(kp, srt=45, cex=2)
  }
  
  if(chr != "canonical"){
    kp <- plotKaryotype(plot.type = 1,   plot.params = pp, chromosomes = chr, cex= 1.5, genome = "hg38")
    kpAddBaseNumbers(kp, 20e6, cex = 2)
    point.cx <- 0.2
    
  }
  kpPlotMarkers(kp, data = gene.markers.gf, labels = gene.markers.gf$Genes, cex =1.70, lwd = 1, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1.35, ignore.chromosome.ends = TRUE,
                marker.parts = c(1.2, 0.01, 0.01))
  
  
  for(i in seq_len(length(sample.names))){
    sn <- names(sample.names)[i]
    
    # image.dir <- file.path(execution.dir,"results",sample.names[i],"images")
    
    #sample tracks
    at <- autotrack(current.track = i, total.tracks = 8, margin = 0.20)
    kpAddLabels(kp, labels = sn, side= "left", r0 =at$r0, r1 = at$r1, cex =2, label.margin = 0.025 )
    
    
    #Gap tracks
    ats <- autotrack(current.track = 1:2, total.tracks = 4,r0 = at$r0,r1=at$r1, margin = 0.20)
    plotCopyNumberCalls(kp, segment.files[[sn]], r0 = ats$r0, r1=ats$r1,  
                        labels = "",  loh.height = 0.2, loh.column = "loh", cn.colors = transparent(cn.colors), track.margin = 0.1)
    
    kpAddLabels(kp, labels = "SNParray",side= "right", r0 =ats$r0, r1 = ats$r1, cex = 2, label.margin = 0.01)
    
    ats <- autotrack(current.track = 1, total.tracks = 4,r0 = at$r0,r1=at$r1, margin = 0.20)
    plotLRR(kp, snps =snp.files[[sn]], labels = "", label.srt = 0, axis.cex =1 , points.cex = point.cx, r0 = ats$r0, r1=ats$r1, track.margin = 0.2,tick.len = 0.02)
    
    ats <- autotrack(current.track = 2, total.tracks = 4,r0 = at$r0,r1=at$r1, margin = 0.20)
    plotBAF(kp,snps =snp.files[[sn]], baf.column = "baf",labels = "",  label.srt = 0, axis.cex =1, points.cex = point.cx, r0 = ats$r0, r1=ats$r1, track.margin = 0.2,tick.len = 0.02 )
    
    #CNVkit tracks
    ats <- autotrack(current.track = 3:4, total.tracks = 4,r0 = at$r0,r1=at$r1, margin = 0.20)
    plotCopyNumberCalls(kp, sample.cn[[sn]], r0 = ats$r0, r1=ats$r1, label.cex = 2.3, 
                        labels = "", label.margin = 0.03,  loh.height = 0, loh.column = "loh", cn.colors = transparent(cn.colors))
    kpAddLabels(kp, labels = "WGS", side= "right", r0 =ats$r0, r1 = ats$r1, cex = 2, label.margin = 0.01 )
    
    ats <- autotrack(current.track = 3, total.tracks = 4, margin = 0.20, r0=at$r0, r1=at$r1)
    plotLRR(kp, snps =sample.log[[sn]], labels = "", axis.cex = 1 , points.cex = point.cx, r0 = ats$r0, r1=ats$r1, line.at.0 = FALSE, track.margin = 0.2,tick.len = 0.02)
    plotCopyNumberCallsAsLines(kp, sample.cn[[sn]], cn.column = "segment.value", style = "segments", col = "red", r0 = ats$r0, r1=ats$r1, ymin = -4, ymax=2, add.axis = FALSE, labels = NA, track.margin = 0.2)
    
    ats <- autotrack(current.track = 4, total.tracks = 4, r0=at$r0, r1=at$r1, margin = 0.20)
    plotBAF(kp,snps =mpnst.snp.lists[[sample.names[i]]], baf.column = "baf",labels = "", axis.cex =1, points.cex = point.cx, r0 = ats$r0, r1=ats$r1, track.margin = 0.2,tick.len = 0.02 )
  }

  dev.off()
  
}


##################### Inversions and translocations of specifc chr ##############################
#SVs of specific genes 
# complete.table.gr <- toGRanges(df.complete)
# 
# chrs <- c(unique(as.character(seqnames(complete.table.gr))))
# 
# # Plot parameters
# cn.colors <- getCopyNumberColors()
# cn.colors[3] <- "white"
# plotDefaultPlotParams(plot.type = 1,plot.params = pp)
# 
# j=1
# for(j in seq_len(length(sample.names)) ){
#   sample <- sample.names[j]
#   ts <- complete.table.gr[complete.table.gr$Sample == names(sample)]
#   transloc <- ts[ts$Variant.type %in%"ICR"]
#   inv <- ts[ts$Variant.type %in%c("INV")]
#   
#   if(!file.exists(file.path(execution.dir,"results",sample.names[j],"Lumpy/Images"))){
#     dir.create(file.path(execution.dir,"results",sample.names[j],"Lumpy/Images"))
#   }
#   
#   for (chr in chrs){
#     # SVs of a specific chr
#     
#     trans.end <-c()
#     transloc.chr <-transloc[seqnames(transloc) == chr]
#     if(length(transloc.chr)>0){
#       trans.end <-toGRanges(transloc.chr$End.bkpoint)
#       
#     }
#     
#     inv.end <-c()
#     inv.chr <-inv[seqnames(inv) == chr]
#     if(length(inv.chr)>0){
#       inv.end <-toGRanges(inv.chr$End.bkpoint)
#     }
#     
#     if(is.null(inv.end) && is.null(trans.end))next
#     
#     png(filename=file.path(execution.dir,"results",sample.names[j],"Lumpy/Images", paste0(sample, "_", chr,"_Translocations&Inv_specific_Chrs_BAF_CN_markers_inv_FILTERED_mw6_test.2.png")), width = 2500, height = 800)
#     
#     pp <- getDefaultPlotParams(plot.type = 2)
#     pp$topmargin <-50
#     pp$bottommargin <- 20
#     pp$data1inmargin <- 30
#     pp$ideogramheight <- 15
#     pp$leftmargin <- 0.15
#     pp$rightmargin <- 0.05
#     pp$data1outmargin <- 100
#     pp$data2height <-20
#     
#     
#     if(length(trans.end)>0){
#       kp <- plotKaryotype(genome = genome, plot.type = 2, chromosomes = c(chr,seqlevels(trans.end)) , cex= 5, plot.params = pp)
#       kpAddBaseNumbers(kp, cex = 2)
#       axis.cex = 2
#     }else if(length(inv.end)>0){
#       kp <- plotKaryotype(genome = genome, plot.type = 2, chromosomes = c(chr) , cex= 5, plot.params = pp)
#       kpAddBaseNumbers(kp, cex = 2)
#       axis.cex = 2
#       
#     }
#     
#     
#     plotCopyNumberCalls(kp, cn.calls = sample.cn[[names(sample)]], data.panel = 1, cn.colors = transparent(cn.colors),  loh.height = 0, loh.column = "loh", labels = NA, label.cex = 5)
#     
#     at <- autotrack(current.track = 2,total.tracks = 2,margin = 0.2)
#     plotBAF(kp,snps = mpnst.snp.lists[[sample.names[j]]], tick.len =0.01, data.panel = 1, axis.cex = axis.cex, label.cex =3.5,  points.cex = 0.5, r0=at$r0, r1=at$r1,labels = "" )
#     at <- autotrack(current.track = 1,total.tracks = 2,margin = 0.2)
#     plotLRR(kp, sample.log[[names(sample)]], lrr.column = "lrr", tick.len =0.01, axis.cex = axis.cex, points.cex = 0.5, r0=at$r0, r1=at$r1, labels = "")
#     
#     
#     # kpPlotMarkers(kp, data = gene.markers.gf, text.orientation = "horizontal", labels = gene.markers.gf$Genes, cex =2, lwd = 2, label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.01, 0.001),ignore.chromosome.ends = TRUE,clipping = FALSE)
#     kpPlotMarkers(kp, data = gene.markers.gf, labels =gene.markers.gf$Genes, cex =2, lwd = 2, srt = 45, label.dist = 0.005, max.iter = 1000, r0= 0, r1 = 1.4,ignore.chromosome.ends = TRUE,
#                   marker.parts = c(0.9, 0.05, 0.05))
#     
#     if(length(transloc.chr) >0){
#       kpPlotLinks(kp,transloc.chr,trans.end, lwd = 4, col = "red")
#       
#     }
#     
#     if(length(inv.chr) >0){      
#       kpPlotLinks(kp,inv.chr, inv.end, lwd = 4, col = "darkblue",r0 = 0, r1=1.33)
#       
#     }
#     dev.off()
#     
#   }
#   
#   
# }
# 
# # Representing All SV of a sample 
# pp <- getDefaultPlotParams(plot.type = 1)
# pp$topmargin <- 30
# pp$bottommargin <- 20
# pp$data1inmargin <- 20
# pp$ideogramheight <- 10
# pp$leftmargin <- 0.15
# pp$rightmargin <- 0.05
# pp$data1outmargin <- 30
# chrs <- c("canonical", unique(as.character(seqnames(sample.cn$`90-8-TL`))))
# j = 2
# chr <- "chr9"
# for(j in seq_len(length(sample.names)) ){
#   sample <- sample.names[j]
#   ts <- filt.list[[sample]][["AllRegionsFiltered"]]
#   transloc <- ts[ts$SVTYPE %in%"ICR"]
#   inv <- ts[ts$SVTYPE %in%c("INV")]
#   
#   if(!file.exists(file.path(execution.dir,"results",sample.names[j],"Lumpy/Images"))){
#     dir.create(file.path(execution.dir,"results",sample.names[j],"Lumpy/Images"))
#   }
#   
#   for (chr in chrs){
#     # chr="canonical"
#     if(chr =="canonical"){
#       
#       ##### All SV of a sample ####
#       # trans.end <- toGRanges(all.sv.samples[[sample]][["ICR"]]$END)
#       trans.end <-c()
#       if(length(transloc)>0){
#         # trans.end <- toGRanges(transloc[transloc$Sample ==sample.names[j]]$END)
#         trans.end <- toGRanges(transloc$END)
#         
#       }
#       
#       inv.end <-c()
#       if(length(inv)>0){
#         # inv.end <- toGRanges(all.sv.samples[[sample]][["INV"]]$END)
#         # inv.end <- toGRanges(inv[inv$Sample ==sample.names[j]]$END)
#         inv.end <- toGRanges(inv$END)
#       }
#       
#       
#     }else{
#       # SVs of a specific chr
#       
#       # transloc.chr <- all.sv.samples[[sample]][["ICR"]][seqnames(all.sv.samples[[sample]][["ICR"]]) == chr]
#       # transloc.chr <-transloc[transloc$Sample==sample.names[j]][seqnames(transloc[transloc$Sample==sample.names[j]]) == chr]
#       # transloc.chr <-transoc.s462[seqnames(transoc.s462) == chr]
#       trans.end <-c()
#       transloc.chr <-transloc[seqnames(transloc) == chr]
#       if(length(transloc.chr)>0){
#         trans.end <-toGRanges(transloc.chr$END)
#         
#       }
#       
#       
#       
#       # inv.chr <- all.sv.samples[[sample]][["INV"]][seqnames(all.sv.samples[[sample]][["INV"]]) == chr]
#       # inv.chr <-inv[inv$Sample==sample.names[j]][seqnames(inv[inv$Sample==sample.names[j]]) == chr]
#       inv.end <-c()
#       inv.chr <-inv[seqnames(inv) == chr]
#       if(length(inv.chr)>0){
#         inv.end <-toGRanges(inv.chr$END)
#       }
#       
#     }
#     
#     if(is.null(inv.end) && is.null(trans.end))next
#     
#     
#     
#     # png(filename=file.path(execution.dir,"results",sample.names[j],"Lumpy/Images", paste0(sample, "_", chr,"_Translocations&Inv_specific_Chrs_BAF_CN_markers_inv_FILTERED_mw6_CTF.png")), width = 2500, height = 1000)
#     
#     if(chr == "canonical"){
#       png(filename=file.path(execution.dir,"results",sample.names[j],"Lumpy/Images", paste0(sample, "_", chr,"_AllTranslocations&Inv_specific_Chrs_BAF_CN_markers_inv_FILTERED_mw6_2.png")), width = 4500, height = 3500)
#       
#       kp <- plotKaryotype(genome = genome, plot.type = 1, chromosomes = chr , plot.params = pp, cex= 3, main = sample)
#       axis.cex = 1
#       
#     }else{
#       png(filename=file.path(execution.dir,"results",sample.names[j],"Lumpy/Images", paste0(sample, "_", chr,"_AllTranslocations&Inv_specific_Chrs_BAF_CN_markers_inv_FILTERED_mw6_2.png")), width = 2500, height = 1500)
#       
#       pp <- getDefaultPlotParams(plot.type = 2)
#       pp$topmargin <-50
#       pp$bottommargin <- 5
#       pp$data1inmargin <- 20
#       pp$ideogramheight <- 15
#       pp$leftmargin <- 0.15
#       pp$rightmargin <- 0.05
#       pp$data1outmargin <- 20
#       pp$data2height <-0
#       pp$data1outmargin <- 100
#       
#       if(length(trans.end)>0){
#         kp <- plotKaryotype(genome = genome, plot.type = 2, chromosomes = c(chr,seqlevels(trans.end)) , cex= 5, plot.params = pp)
#         kpAddBaseNumbers(kp, cex = 2)
#         axis.cex = 1.5
#       }else if(length(inv.end)>0){
#         kp <- plotKaryotype(genome = genome, plot.type = 2, chromosomes = c(chr) , cex= 5, plot.params = pp)
#         kpAddBaseNumbers(kp, cex = 2)
#         axis.cex = 1.5
#         
#       }
#     }
#     
#     # kp <- plotKaryotype(genome = "hg38",chromosomes = chr, plot.params = pp, plot.type = 3, cex = 1.5 )
#     # kpAddCytobandsAsLine(kp,lwd = 2)
#     kpPlotMarkers(kp, data = gene.markers.gf, labels = gene.markers.gf$Genes, cex =2, lwd = 2, srt = 45, label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.01, 0.001),ignore.chromosome.ends = TRUE,clipping = FALSE)
#     
#     plotCopyNumberCalls(kp, cn.calls = sample.cn[[names(sample)]], data.panel = 1, cn.colors = transparent(cn.colors),  loh.height = 0, loh.column = "loh", labels = NA, label.cex = 5)
#     
#     at <- autotrack(current.track = 2,total.tracks = 2,margin = 0.2)
#     plotBAF(kp,snps = mpnst.snp.lists[[sample.names[j]]], tick.len =0.01, data.panel = 1, axis.cex = axis.cex, label.cex =3.5,  points.cex = 0.5, r0=at$r0, r1=at$r1,labels = "" )
#     at <- autotrack(current.track = 1,total.tracks = 2,margin = 0.2)
#     plotLRR(kp, sample.log[[names(sample)]], lrr.column = "lrr", tick.len =0.01, axis.cex = axis.cex, points.cex = 0.5, r0=at$r0, r1=at$r1, labels = "")
#     
#     if(chr == "canonical"){
#       
#       if(length(transloc)>0){
#         kpPlotLinks(kp, transloc,trans.end, lwd = 4, col = "red")
#         
#       }
#       
#       if(length(inv)>0){
#         kpPlotLinks(kp, inv, inv.end, lwd = 4, col = "darkblue",r0 = 0, r1=1.33)
#         
#       }
#       
#     }else{
#       
#       if(length(transloc.chr) >0){
#         kpPlotLinks(kp,transloc.chr,trans.end, lwd = 4, col = "red")
#         
#       }
#       
#       if(length(inv.chr) >0){      
#         kpPlotLinks(kp,inv.chr, inv.end, lwd = 4, col = "darkblue",r0 = 0, r1=1.33)
#         
#       }
#       
#     }
#     
#     dev.off()
#     
#   }
#   
# }

#########################################################
#             Chromosomal representation                #
#########################################################
# # Plot parameters
# cn.colors <- getCopyNumberColors()
# cn.colors[3] <- "white"
# 
# pp <- getDefaultPlotParams(plot.type = 1)
# pp$topmargin <- 30
# pp$bottommargin <- 20
# pp$data1inmargin <- 5
# pp$ideogramheight <- 2
# pp$ideogramlateralmargin <- 0.008
# pp$leftmargin <- 0.1
# # sample <- "S462"
# chrs <- c("canonical",unique(as.character(seqnames(sample.cn[["S462"]]))))
# chr = "canonical"
# sample = "NMS.2"
# i=3
# for(i in seq_len(length(sample.names)) ){
#   sample <- sample.names[i]
#   for (chr in chrs){
#     cnvkit.dir <- file.path(execution.dir,"results",sample,"CNVkit/ExcludedRegions2_ploidy/images")
#     if(!file.exists(cnvkit.dir))dir.create(cnvkit.dir)
#     png(filename=file.path(cnvkit.dir, paste0(sample, "_", chr,"_BAF_CN_markers.png")), width = 2500, height = 1000)
#     # png(filename=file.path(execution.dir,"results",sample,"CNVkit/images", paste0(sample, "_", chr,"_BAF_CN_markers.png")), width = 2500, height = 1000)
#     
#     kp <- plotKaryotype(genome = "hg38",chromosomes = chr, plot.params = pp, plot.type = 4, cex = 1.5 )
#     # kpAddCytobandsAsLine(kp,lwd = 2)
#     if(chr!="canonical")kpAddBaseNumbers(kp, cex =1.5)
#     
#     kpPlotMarkers(kp, data = gene.markers.gf, labels = gene.markers.gf$Genes, cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1500, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.01, 0.001),ignore.chromosome.ends = TRUE,clipping = FALSE)
#     
#     at <- autotrack(current.track = 1, total.tracks = 2,margin = 0.20)
#     plotCopyNumberCalls(kp, sample.cn[[names(sample)]], r0 = at$r0, r1=at$r1, label.cex = 2.3, 
#                         labels = NA,  loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
#     
#     plotLRR(kp, snps =sample.log[[names(sample)]],  lrr.column = "lrr",labels = "LRR",label.cex = 2.3, label.srt = 0, label.margin = 0.08, axis.cex =1.8 , points.cex = 0.5, r0 = at$r0, r1=at$r1)
#     plotCopyNumberCallsAsLines(kp, sample.cn[[names(sample)]], cn.column = "segment.value", style = "segments", col = "red", r0 = at$r0, r1=at$r1, ymin = -4, ymax=2, add.axis = FALSE, labels = NA)
#     
#     at <- autotrack(current.track = 2, total.tracks = 2)
#     plotCopyNumberCalls(kp, sample.cn[[names(sample)]], r0 = at$r0, r1=at$r1, label.cex = 2.3,
#                         labels = NA,  loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
#     
#     plotBAF(kp,snps = mpnst.snp.lists[[sample]][[sample]], baf.column = "baf",labels = "BAF", label.cex = 2.3, label.srt = 0, label.margin = 0.08, axis.cex =1.8, points.cex = 0.5, r0 = at$r0, r1=at$r1 )
#     
#     
#     
#     # gene.markers <- c(gene.markers,toGRanges("chr22:20982296-20999031"),toGRanges("chr22:29603555-29698599"),toGRanges("chr22:23786965-23838008"), toGRanges("chr7:148807384-148884248"))
#     # gene.markers$ssymbol[c(6:9)]<- c("LZTR1","NF2","SMARCB1", "EZH2")
#     
#     
#     
#     
#     dev.off()
#     
#   }
# }



