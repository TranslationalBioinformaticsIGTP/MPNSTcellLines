#####################################################
# Script Supp Document ST88-14 and T265 cell lines  #
#####################################################
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
execution.dir <- "./MPNST_cellLines/SNP-array"

# Loading the SNParray data from WGS
load(file.path(execution.dir,"results","MPNST_WGS_cellLines_BAF_Filtered.RData"))

##### loading GAP SNP-array data:download data from synapse ####
data.dir <- ""
result.dir <- "./MPNST_cellLines/SNP-array/results/paperCDKN2A"

cell.lines <- read.table("./MPNST_cellLines/SNP-array/cellLinesDataBernatDir.csv",sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sample.names <- cell.lines$sample.name
# cell.lines <-cell.lines[c(4,2,1,5,3,6,8,7),]
rownames(cell.lines)<- sample.names
sample.names

####Loading segment data
segment.files <- list()
snp.files <- list()
i=13
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


# Plot parameters
###################
####Plot BAF LRR CN images
load(file = file.path("./MPNST_cellLines/WGS","results","table.to.plot.RData"))


cn.colors <- getCopyNumberColors()
cn.colors[3] <- "white"

pp <- getDefaultPlotParams(plot.type = 1)
pp$topmargin <- 30
pp$bottommargin <- 30
pp$data1inmargin <- 5
pp$ideogramheight <- 2
pp$ideogramlateralmargin <- 0.008
pp$leftmargin <- 0.1

chr <- "chr16"
i=13
for(i in seq_len(nrow(cell.lines))){
  cl <- rownames(cell.lines)[i]
  
  result.dir <- file.path("./MPNST_cellLines/SNP-array/results/",cl)
  if(!file.exists(result.dir))dir.create(result.dir)
  
  png(file.path(result.dir, paste0(cl,"_BAF_LRR_CN_",chr,".3.png")), width = 2500, height = 800)
  if(chr == "canonical"){
    kp <- plotKaryotype(plot.type = 4,   plot.params = pp, chromosomes = chr, cex= 1.35, genome = "hg38")
    # kpAddCytobandsAsLine(kp)
    point.cx <- 0.2
    # kpAddBaseNumbers(kp,cex=1)
  }else{
    kp <- plotKaryotype(plot.type = 1,  main = cl, plot.params = pp, chromosomes = chr, cex= 3, genome = "hg38")
    kpAddBaseNumbers(kp, 20e6, cex = 2)
    
  }
  # kp <- plotKaryotype(plot.type = 1,  main = cl, plot.params = pp, chromosomes = chr, cex= 3, genome = "hg38")
  # kpAddBaseNumbers(kp, 20e6, cex = 2)
  kpPlotMarkers(kp, data = table.to.plot[["90-8-TL"]], labels =table.to.plot[["90-8-TL"]]$Genes ,label.color = "black",  cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1000, r0= 0, r1 = 1.4,ignore.chromosome.ends = TRUE,
                marker.parts = c(0.99, 0.03, 0.01))
  
  plotCopyNumberCalls(kp, segment.files[[i]], r0 = 0, r1=0.10, label.cex = 2.3, label.margin = 0.005, 
                      labels = "CN-SNParray", track.margin = 0.1, loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
  # plotCopyNumberCalls(kp, segment.files, r0 = 0, r1=1, label.cex = 2.3, label.margin = 0.005, 
  #                     labels = NULL, track.margin = 0.1, loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
  
  # at <- autotrack(current.track = 1, total.tracks = 2,margin = 0.20)
  plotLRR(kp, snps =snp.files[[rownames(cell.lines)[i]]], ymin = -2, labels = "LRR",label.cex = 2.3, label.srt = 0, label.margin = 0.05, axis.cex =1.6 , points.cex = point.cx, r0 =0.12, r1=0.57, track.margin = 0.25, line.at.0.col = "red")

  # at <- autotrack(current.track = 2, total.tracks = 2)
  # plotCopyNumberCalls(kp, segment.files[[cl]], r0 = at$r0, r1=at$r1, label.cex = 2.3,label.margin = 0.07, 
  # labels = NA, track.margin = 0.1, loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
  # at <- autotrack(current.track = 2, total.tracks = 2)
  plotBAF(kp,snps =snp.files[[rownames(cell.lines)[i]]], baf.column = "baf",labels = "BAF", label.cex = 2.3, label.srt = 0, label.margin = 0.05, axis.cex =1.6, points.cex = point.cx, r0 = 0.59, r1=1, track.margin = 0.25 )
  
  
  # kpPlotMarkers(kp, data =gene.markers, r1 = 1.45, text.orientation = "horizontal", labels = gene.markers$symbol, cex = 3, lwd = 2)
  # kpPlotMarkers(kp, data = gene.markers, labels = gene.markers$symbol, cex =3.2, lwd = 2, text.orientation = "horizontal", label.dist = 0.045, max.iter = 1000, r0= 0, r1 = 1,y = 1.05, marker.parts = c(0.98, 0.4, 0.1))
  
  
  # kpPlotLinks(kp,transloc.chr9,trans.end, lwd = 5, col = "orange")
  # kpPlotLinks(kp,inv.chr9, inv.end, lwd = 8, col = "darkblue",r0 = 0, r1=1.33)
  
  dev.off()
}


### Sup Doc 1
sup.doc.1 <- segment.files[c(2,11,10,12,9,13,14)]
s
names(sup.doc.1)
cell.names <- names(sup.doc.1)
names(cell.names) <- c("ST88-14_Lab1", "ST88-14_PT_JF",
                       "ST94531_JF", "ST88-14_Lab2",
                       "T265_GdV","T265_Lab2","ST88-14_CLp69_JF")

regions.dif  <- toGRanges(c("chr2:68000000-72000000","chr3:1-198295559","chr5:1-181538259",
                            "chr6:1-58000000","chr7:1-159345973","chr8:1-145138636","chr9:68000000-138394717",
                            "chr13:54000000-114364328","chr19:1-58617616"))
cn.colors <- getCopyNumberColors()
cn.colors[3] <- "white"

pp <- getDefaultPlotParams(plot.type = 3)
pp$topmargin <-33
# pp$topmargin <-5 #no genes
pp$bottommargin <- 12
pp$data1inmargin <- 1
pp$ideogramheight <- 2
pp$ideogramlateralmargin <- 0.008
pp$leftmargin <- 0.2
pp$data2height <- 0
pp$data1outmargin <- 2
plotDefaultPlotParams(plot.type = 3, plot.params = pp)

chr <- "canonical"

result.dir <- file.path("./MPNST_cellLines/SNP-array/results/Group_images")
if(!file.exists(result.dir))dir.create(result.dir)
  
png(file.path(result.dir, paste0("SupDoc1.1_BAF_LRR_CN_",chr,".Genes.png")), width = 1800, height = 1000)
kp <- plotKaryotype(plot.type = 4, plot.params = pp, chromosomes = chr, cex= 1, genome = "hg38", labels.plotter = NULL)
kpAddChromosomeNames(kp, srt=45, cex=2,yoffset = 0)
kpPlotMarkers(kp, data = table.to.plot[["90-8-TL"]], labels =table.to.plot[["90-8-TL"]]$Genes,line.color = "#CCCCCC",label.color = "black",  cex =1.8, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1000, r0= 0, r1 = 1.4,ignore.chromosome.ends = TRUE,
              marker.parts = c(0.99, 0.03, 0.01))
#Diff (88-14 Ref vs JF)
# kpPlotRegions(kp,data = regions.dif,col = transparent("#eed7f4ff"),border = transparent("#eed7f4ff"))
for(i in seq_len(length(cell.names[c(1:3,7)]))){
# for(i in seq_len(length(cell.names[c(1,5,6)]))){
# for(i in seq_len(length(cell.names[c(1,4,5,6,7)]))){
  # cl <- names(cell.names)[c(1,5,6)][i]
  # sn <- cell.names[c(1,5,6)][i]
  cl <- names(cell.names)[c(1:3,7)][i]
  sn <- cell.names[c(1:3,7)][i]
  # cl <- names(cell.names)[c(1,4,5,6,7)][i]
  # sn <- cell.names[c(1,4,5,6,7)][i]
  point.cx <- 0.2

  #F1 supDoc1
  if(cl =="ST88-14_Lab1"){
    r0=0.775
    r1 =1
  }else if(cl =="ST88-14_CLp69_JF"){
    r0=0.452
    r1=0.677
  }else if(cl == "ST88-14_PT_JF"){
    r0=0.227
    r1 =0.450

  }else if(cl == "ST94531_JF"){
    r0=0
    r1 =0.225
  }

  # #F2 SupDoc1
  # if(cl =="ST88-14_Lab1"){
  #   r0=0.7
  #   r1 =1
  # }else if(cl == "T265_GdV"){
  #   r0=0.32
  #   r1 =0.62
  # }else if(cl == "T265_Lab2"){
  #   r0=0
  #   r1 =0.3
  # }
  
  # # F2 supDoc2
  # if(cl =="ST88-14_Lab1"){
  #   r0=0.82
  #   r1 =1
  # }else if(cl =="ST88-14_CLp69_JF"){
  #   r0=0.56
  #   r1=0.74
  # }else if(cl == "ST88-14_Lab2"){
  #   r0=0.36
  #   r1 =0.54
  # }else if(cl == "T265_Lab2"){
  #   r0=0.18
  #   r1 =0.34
  # }else if(cl == "T265_GdV"){
  #   r0=0
  #   r1 =0.16
  # }
  # at <- autotrack()
    plotCopyNumberCalls(kp, segment.files[[sn]], r0 = r0, r1=r0+0.03, label.cex = 2.3, label.margin = 0.005, 
                        labels = "", track.margin = 0.1, loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
    
    kpAddLabels(kp, labels = cl, side= "left", r0 =r0, r1 = r1+0.03, cex =2.3, label.margin = 0.060,pos=2)
    
    ats <- autotrack(current.track = 1, total.tracks = 2, r0 = r0+0.03, r1=r1, margin = 0.2)
    
    plotLRR(kp, snps =snp.files[[sn]], ymin = -2, r0 = ats$r0, r1=ats$r1, labels = "LRR",label.cex = 1.5, label.srt = 90, label.margin = 0.025, axis.cex =1.6 , points.cex = point.cx, track.margin = 0.25, line.at.0.col = "red")
    ats <- autotrack(current.track = 2, total.tracks = 2,r0 = r0+0.03, r1=r1, margin = 0.2)
    
    plotBAF(kp,snps =snp.files[[sn]], baf.column = "baf", r0 = ats$r0, r1=ats$r1, labels = "BAF", label.cex = 1.5, label.srt = 90, label.margin = 0.025, axis.cex =1.6, points.cex = point.cx, track.margin = 0.25 )
    
  
}
  dev.off()

  
  ### Sup Doc 2
  sup.doc.2 <- segment.files[c(4,15,16)]
  s
  names(segment.files)
  cell.names <- names(sup.doc.2)
  names(cell.names) <- c("S462_Lab3", "S462_Lab2",
                         "S462_Lab4")
  
  regions.dif <-toGRanges(c("chr1:1-44000000","chr4:55000000-190214555","chr5:1-135000000","chr7:1-159345973",
                            "chr9:1-138394717","chr16:1-35000000"))
  # toDataframe(kp$genome)
  cn.colors <- getCopyNumberColors()
  cn.colors[3] <- "white"
  
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$topmargin <-33
  # pp$topmargin <-5 #no genes
  pp$bottommargin <- 12
  pp$data1inmargin <- 1
  pp$ideogramheight <- 2
  pp$ideogramlateralmargin <- 0.008
  pp$leftmargin <- 0.2
  pp$data2height <- 0
  pp$data1outmargin <- 2
  plotDefaultPlotParams(plot.type = 3, plot.params = pp)
  
  chr <- "canonical"
  
  result.dir <- file.path("./MPNST_cellLines/SNP-array/results/Group_images")
  if(!file.exists(result.dir))dir.create(result.dir)
  
  png(file.path(result.dir, paste0("SupDoc2.1_BAF_LRR_CN_",chr,".Genes.shadow.png")), width = 1800, height = 1000)
  kp <- plotKaryotype(plot.type = 4, plot.params = pp, chromosomes = chr, cex= 1, genome = "hg38", labels.plotter = NULL)
  kpAddChromosomeNames(kp, srt=45, cex=2,yoffset = 0)
  kpPlotMarkers(kp, data = table.to.plot[["90-8-TL"]], labels =table.to.plot[["90-8-TL"]]$Genes ,line.color = "#CCCCCC",label.color = "black",  cex =1.8, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1000, r0= 0, r1 = 1.4,ignore.chromosome.ends = TRUE,
                marker.parts = c(0.99, 0.03, 0.01))
  kpPlotRegions(kp,data = regions.dif,col = transparent("#eed7f4ff"),border = transparent("#eed7f4ff"))
  
  for(i in seq_len(length(cell.names))){
    
    cl <- names(cell.names)[i]
    sn <- cell.names[i]
    
    point.cx <- 0.2
 
    # F1 SupDoc2
    if(cl =="S462_Lab3"){
      r0=0.7
      r1 =1
    }else if(cl == "S462_Lab2"){
      r0=0.32
      r1 =0.62
    }else if(cl == "S462_Lab4"){
      r0=0
      r1 =0.3
    }
    
    plotCopyNumberCalls(kp, segment.files[[sn]], r0 = r0, r1=r0+0.03, label.cex = 2.3, label.margin = 0.005, 
                        labels = "", track.margin = 0.1, loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
    
    kpAddLabels(kp, labels = cl, side= "left", r0 =r0, r1 = r1+0.03, cex =2.3, label.margin = 0.060,pos=2)
    
    ats <- autotrack(current.track = 1, total.tracks = 2, r0 = r0+0.03, r1=r1, margin = 0.2)
    
    plotLRR(kp, snps =snp.files[[sn]], ymin = -2, r0 = ats$r0, r1=ats$r1, labels = "LRR",label.cex = 1.5, label.srt = 90, label.margin = 0.025, axis.cex =1.6 , points.cex = point.cx, track.margin = 0.25, line.at.0.col = "red")
    ats <- autotrack(current.track = 2, total.tracks = 2,r0 = r0+0.03, r1=r1, margin = 0.2)
    
    plotBAF(kp,snps =snp.files[[sn]], baf.column = "baf", r0 = ats$r0, r1=ats$r1, labels = "BAF", label.cex = 1.5, label.srt = 90, label.margin = 0.025, axis.cex =1.6, points.cex = point.cx, track.margin = 0.25 )
    
    
  }
  dev.off()
  
  ################################################################################################################
  #                                           SnpArray: MPNSTs CellLines                                        #
  ################################################################################################################
  
  
  # #######################Plot BAF LRR CN images
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
  # 
  # chr <- "canonical"
  # i=1
  # for(i in seq_len(nrow(cell.lines))){
  #   cl <- rownames(cell.lines)[i]
  #   result.dir <- file.path("./MPNST_cellLines/SNP-array/results/",cl)
  #   if(!file.exists(result.dir))dir.create(result.dir)
  #   
  #   png(file.path(result.dir, paste0(name.cells[i],"_BAF_LRR_CN_",chr,".3.png")), width = 2500, height = 800)
  #   if(chr == "canonical"){
  #     kp <- plotKaryotype(plot.type = 4,   plot.params = pp, chromosomes = chr, cex= 1.35, genome = "hg38")
  #     # kpAddCytobandsAsLine(kp)
  #     point.cx <- 0.2
  #     # kpAddBaseNumbers(kp,cex=1)
  #   }else{
  #     kp <- plotKaryotype(plot.type = 1,  main = name.cells[i], plot.params = pp, chromosomes = chr, cex= 3, genome = "hg38")
  #     kpAddBaseNumbers(kp, 20e6, cex = 2)
  #     
  #   }
  #     kpPlotMarkers(kp, data = table.to.plot[[sn]], labels =table.to.plot[[sn]]$Genes ,label.color = "black",  cex =2, lwd = 2, text.orientation = "vertical", label.dist = 0.005, max.iter = 1000, r0= 0, r1 = 1.4,ignore.chromosome.ends = TRUE,
  #                 marker.parts = c(0.99, 0.03, 0.01))
  #   
  #   plotCopyNumberCalls(kp, segment.files[[i]], r0 = 0, r1=0.10, label.cex = 2.3, label.margin = 0.005, 
  #                       labels = "CN-SNParray", track.margin = 0.1, loh.height = 0.2, loh.column = "loh", cn.colors = cn.colors)
  #     plotLRR(kp, snps =snp.files[[cl]], ymin = -2, labels = "LRR",label.cex = 2.3, label.srt = 0, label.margin = 0.05, axis.cex =1.6 , points.cex = point.cx, r0 =0.12, r1=0.57, track.margin = 0.25)
  #   
  #      plotBAF(kp,snps =snp.files[[cl]], baf.column = "baf",labels = "BAF", label.cex = 2.3, label.srt = 0, label.margin = 0.05, axis.cex =1.6, points.cex = point.cx, r0 = 0.59, r1=1, track.margin = 0.25 )
  # 
  #   
  #   dev.off()
  # }
  # 
  # 