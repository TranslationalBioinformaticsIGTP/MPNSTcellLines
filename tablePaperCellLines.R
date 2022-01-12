###########################################################
# Script to generate a table summaryzing the alterations ## 
#########  on specific genes of the cell lines  ###########
###########################################################

# Packages needed
library(karyoploteR)
library(CopyNumberPlots)

#Directories 
execution.dir <- "./MPNST_cellLines/WGS/"
results.dir <- file.path(execution.dir, "results")
table.dir <- file.path(results.dir, "tableSNVsSVs_PaperCellLines")
if(!file.exists(table.dir))dir.create(table.dir)
"chr7:148,807,385-148,884,245"

# Loading data
# load(file = file.path(results.dir,"table.to.plot.RData"))
load(file = file.path(results.dir,"table.to.plot.toKeep.3.RData"))

table.plot <- toDataframe(unlist(as(table.to.plot,"GRangesList")))
cells <- unique(table.plot$Sample)
cells <- c(cells,"sp")
cells <- cells[c(1:5,9,6:8)]
# genes <- c("NF1", "CDKN2A", "SUZ12", "EED", "sp", "TP53", "PTEN", "RB", "sp", "BIRC5", "G1", "G2", "G3", "sp", "TWIST1", "G4", "G5", "G6")
genes <- c("NF1","CDKN2A","sp", "SUZ12", "EED", "sp",
  "TP53", "PTEN","RB1","sp",
  "EGFR", "ERBB2", "ERBB3", "PDGFRA", "KIT", "MET", "HGF", "CXCR4","sp",
  "BIRC5",  "TWIST1", "SOX9", "AURKA","AURKB")

# genes <- c("NF1","CDKN2A","sp", "SUZ12", "EED", "sp","TP53", "PTEN","RB1")
# genes <- c("EGFR", "ERBB2", "ERBB3", "PDGFRA", "KIT", "MET", "HGF", "CXCR4")
# genes <- c("BIRC5",  "TWIST1", "SOX9", "AURKA","AURKB")
         



# cells <- list("S462", "ST88-14", "90-8", "sp", "STS-26T", "HS-PSS", "HS-Sch2")
# table.alterations <- list()

table.plot.cells <- split(table.plot, table.plot$Sample)

#LOH
loh <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
rownames(loh) <- genes
colnames(loh) <- cells
# loh[2,3] <- TRUE #ejemplo

#SNV
snv <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# snv[1,6] <- TRUE #ejemplo
rownames(snv) <- genes
colnames(snv) <- cells

#SV
sv <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# sv[1,4] <- TRUE #ejemplo
rownames(sv) <- genes
colnames(sv) <- cells

#Hom.loss
hom.loss <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# hom.loss[1,7] <- TRUE #ejemplo
rownames(hom.loss) <- genes
colnames(hom.loss) <- cells

#Gain 
gain <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# gain[4,6] <- TRUE #ejemplo
rownames(gain) <- genes
colnames(gain) <- cells

#het.loss
het.loss <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# het.loss[5,2] <- TRUE #ejemplo
rownames(het.loss) <- genes
colnames(het.loss) <- cells

#inact
inact <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# inact[2,4] <- TRUE #ejemplo
rownames(inact) <- genes
colnames(inact) <- cells

i=1
j=2
for(i in seq_len(length(table.plot.cells))){
  cl <- names(table.plot.cells)[i]
  tb <- table.plot.cells[[cl]]
  tb <-tb[tb$Genes%in% genes,]
  for(j in seq_len(nrow(tb))){
    gene <- tb$Genes[j]
    loh[gene,cl] <- as.logical(tb$loh[j])
    if(is.na(loh[gene,cl])){
      loh[gene,cl] <- FALSE
    }
    #SV and SNV
    if(is.na(tb$variant[j])){
      snv[gene,cl] <- FALSE
      sv[gene,cl] <- FALSE
    }else{
      if(tb$variant[j]=="SNV"){
        snv[gene,cl] <- TRUE
      }
      if(tb$variant[j] == "SV"){
        sv[gene,cl] <- TRUE
      }
    }
    
    
    
    #CNVs
    if(is.na(tb$cn[j])){
      gain[gene,cl] <- FALSE
      hom.loss[gene,cl] <- FALSE
      het.loss[gene,cl] <- FALSE
    }else{
      if(tb$cn[j]=="Gain"){
        gain[gene,cl] <- TRUE
      }
      if(tb$cn[j] == "HomLoss"){
        hom.loss[gene,cl] <- TRUE
      }
      if(tb$cn[j]=="HetLoss"){
        het.loss[gene,cl]<- TRUE
      }
      
    }
    
    
    
    #inactivation
    if(tb$altered[j] == "complete"){
      inact[gene,cl] <- TRUE
    }else{
      inact[gene,cl] <- FALSE
    }
    
  }
}





# 
# 
# 
# #Para los datos lo más fácil es crear una matriz de length(genes) filas y length(cells) columnas para cada cosa
# 
# unique(table.plot$Genes)
# 
# # creating the table.matrix to plot
# 
# 
# 
# 
# 
#Create table

plot(x=20:80, y=20:80)

createCell <- function(x, y, size=25, loh=FALSE, snv=FALSE, snv.col="black", snv.cex=3, sv=FALSE, sv.col="black", sv.cex=2.5, inact=FALSE, gain=FALSE, het.loss=FALSE, hom.loss=FALSE) {
  mid <- ((size-1)/2)

  rect(xleft=x-mid, ybottom = y-mid, xright = x+mid, ytop = y+mid, col="#CCCCCC", border=NA)

  if(loh) {
    rect(xleft=x-mid, xright=x+mid, ybottom=y-mid, ytop=y-(mid/2), col="dodgerblue", border=NA)
  }
  if(snv) {
    points(x, y, cex=snv.cex, pch=16, col=snv.col)
  }
  if(sv) {
    points(x, y, cex=sv.cex, pch=17, col=sv.col)
  }
  if(gain) {
    # rect(xleft=x-mid, ybottom = y-mid, xright = x+mid, ytop = y+mid, col=NA, border="red", lwd=6)
    rect(xleft=x-mid, ybottom = y-mid, xright = x+mid, ytop = y+mid, col=NA, border="#C04343", lwd=3)
    
  }
  if(het.loss) {
    rect(xleft=x-mid, ybottom = y-mid, xright = x+mid, ytop = y+mid, col=NA, border="lightgreen", lwd=3)
  }
  if(hom.loss) {
    rect(xleft=x-mid, ybottom = y-mid, xright = x+mid, ytop = y+mid, col=NA, border="darkgreen", lwd=3)
  }
  if(inact) {
    segments(x0 = x-mid, y0 = y-mid, x1 = x+mid, y1 = y+mid, col="#333333", lwd=6)
    segments(x0=x-mid, y1 = y-mid, x1 = x+mid, y0 = y+mid, col="#333333", lwd=6)

  }
}

# plot(x=20:80, y=20:80, type="n")
# createCell(50, 50, loh=TRUE, sv=TRUE, snv=TRUE, hom.loss=TRUE, inact=TRUE)
# # createCell(50, 50, loh=loh, sv=sv, snv=snv, hom.loss=hom.loss, inact=inact)
# 
# # 
# # genes <- c("NF1", "CDKN2A", "SUZ12", "EED", "sp", "TP53", "PTEN", "RB", "sp", "BIRC5", "G1", "G2", "G3", "sp", "TWIST1", "G4", "G5", "G6")
# # cells <- list("S462", "ST88-14", "90-8", "sp", "STS-26T", "HS-PSS", "HS-Sch2")
# 
# 
# #Para los datos lo más fácil es crear una matriz de length(genes) filas y length(cells) columnas para cada cosa
# 
# #LOH
# loh <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# loh[2,3] <- TRUE #ejemplo
# 
# 
# #SNV
# snv <- matrix(rep(FALSE, length(genes)*length(cells)), nrow=length(genes))
# snv[1,6] <- TRUE #ejemplo
# 





#Create an empty plot

graphics::par(mar=c(0,0,0,0)+0.1)

left.mar <- 100
right.mar <-  40
top.mar <- 100
bottom.mar <- 40
wide.mar <- 10
narrow.mar <- 5

cell.size <- 25

height.total <- top.mar + bottom.mar + length(genes)*(narrow.mar + cell.size) + length(which(genes=="sp"))
width.total <- right.mar + left.mar +  length(cells)*(narrow.mar + cell.size) +  length(which(cells=="sp"))


#png("Test.png", width=2*width.total, height = 2*height.total)
# svg(file.path(execution.dir,"results/tableSNVsSVs_PaperCellLines","Others.svg"), width=width.total/50, height = height.total/50)

#positions of the gene names
gy <- height.total - top.mar + 20
gx <- left.mar - 20

graphics::plot(0, type="n", xlim=c(0,width.total), ylim=c(0, height.total), axes=FALSE, ylab="", xlab="", xaxs="i", yaxs="i")

abline(h=0, v=0)

abline(h=gy, v=gx)

first.g <- height.total - top.mar - cell.size/2
first.c <- left.mar + cell.size/2

cy <- first.g
for(ng in 1:length(genes)) {
  if(genes[ng]=="sp") {
    cy <- cy - wide.mar - narrow.mar
  } else {
    #plot the gene name
    text(x=gx, y=cy, label=genes[ng], pos=2)    
    #plot the cells
    cx <- first.c    
    for(nc in 1:length(cells)) {
      if(cells[nc]=="sp") {
        cx <- cx + wide.mar + narrow.mar
      } else {
        if(ng==1) {
          #plot the cell line name
          text(x=cx, y=gy+25, label=cells[nc], srt=90)
        }
        createCell(x=cx, y=cy, 
                   loh=loh[ng, nc], 
                   snv=snv[ng,nc], 
                   sv = sv[ng,nc],
                   inact = inact[ng,nc],
                   gain = gain[ng,nc],
                   het.loss = het.loss[ng,nc],
                   hom.loss = hom.loss[ng,nc]     )
        #text(x=cx, y=cy, label=paste0(ng, ",", nc))
        cx <- cx + cell.size + narrow.mar
      }
    }
    cy <- cy - cell.size - narrow.mar
  }
}



dev.off()
