##################### Mut Signatures #############
library(mutSignatures)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(BSgenome.Hsapiens.UCSC.hg38)

#Get 30 different colors
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


execution.dir<- "./MPNST_cellLines/WGS"
image.dir <- file.path(execution.dir,"results/MutSignatures/images")
if(!file.exists(image.dir))dir.create(image.dir)
mut.signatures.df <- read.table(file = file.path(execution.dir,"results",
                                                 paste0("AllSamples_MutSignVariants.3.csv")),sep = "\t", header=T,stringsAsFactors = F)
# mut.signatures.df <- mut.signatures.df[mut.signatures.df$SAMLEID != "STS-26T",]
# 
# mut.signatures.df <- read.table(file = file.path(execution.dir,"results", 
#                                                  paste0("AllSamples_MutSignVariants.TMB.csv")),sep = "\t", header=T,stringsAsFactors = F)

x <- filterSNV(dataSet = mut.signatures.df, seq_colNames = c("REF","ALT"))
genome <- BSgenome.Hsapiens.UCSC.hg38

x <- attachContext(mutData = x,
                   chr_colName = "CHROM",
                   start_colName = "POS",
                   end_colName = "POS",
                   nucl_contextN = 3,
                   BSGenomeDb = genome)
# Remove mismatches
x <- removeMismatchMut(mutData = x,                  # input data.frame
                       refMut_colName = "REF",       # column name for ref base
                       context_colName = "context",  # column name for context
                       refMut_format = "N")
# Compute mutType
x <- attachMutType(mutData = x,                      # as above
                   ref_colName = "REF",              # column name for ref base
                   var_colName = "ALT",              # column name for mut base
                   context_colName = "context")

# Mutation Counts
cells.signature.counts <- countMutTypes(mutTable = x,
                             mutType_colName = "mutType",
                             sample_colName = "SAMLEID")

# save(cells.signature.counts, file = file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts.2.RData"))
# load(file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts.2.RData"))
# save(cells.signature.counts, file = file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts_NO_26T.RData"))


df<- mutSignatures::as.data.frame(cells.signature.counts)
colSums(df)
#Get cosmic signatures 
cosmic <-mutSignatures::getCosmicSignatures()
# # 
# # cosmic.selected <- c("COSMIC.16","COSMIC.7","COSMIC.30","COSMIC.12","COSMIC.19","COSMIC.5")
# # cosmic.selected <- c("COSMIC.1","COSMIC.3","COSMIC.5","COSMIC.8","COSMIC.9",
# #                      "COSMIC.12","COSMIC.16","COSMIC.19","COSMIC.20",
# #                      "COSMIC.25","COSMIC.26","COSMIC.30")
# # cosmic.selected <- c("COSMIC.16","COSMIC.7","COSMIC.5")
# # 
# # # cosmic.selected <- rownames(msig1$distanceMatrix)[rowSums(msig1$distanceMatrix)<0.6]
# # cosmic.df <- mutSignatures::as.data.frame(cosmic)
# # cosmic.df <- cosmic.df[,cosmic.selected]
# # cosmic <- mutSignatures::as.mutation.signatures(cosmic.df)
# # 
# 
# # how many signatures should we extract? 
# num.sign <- 2
# load(file.path(execution.dir,"results","MutSignatures",
#                paste0("MPNSTcells_MutSignatureAnalysis_1000Iter_",num.sign,"sign.RData")))
# # 
# # # Define parameters for the non-negative matrix factorization procedure.
# # # you should parallelize if possible
# # params <-
# #   mutSignatures::setMutClusterParams(
# #     num_processesToExtract = num.sign,    # num signatures to extract
# #     num_totIterations = 50,               # bootstrapping: usually 500-1000
# #     num_parallelCores = 3)                # total num of cores to use (parallelization)
# # 
# # # Extract new signatures - may take a while
# # cell.signature.analysis <-
# #   decipherMutationalProcesses(input = cells.signature.counts,
# #                               params = params)
# 
# # # save(cells.signature.counts, file = file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts.RData"))
# # save(cell.signature.analysis, file = file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts_NO_26T.RData"))
# 
# # Retrieve signatures (results)
# cell.signature.sig <- cell.signature.analysis$Results$signatures
# 
# # Retrieve exposures (results)
# cell.signature.exp <- cell.signature.analysis$Results$exposures
# 
# # Plot signature 1 (standard barplot, you can pass extra args such as ylim)
# msigPlot(cell.signature.sig, signature =1 , ylim = c(0, 0.05))
# 
# 
# png(file.path(image.dir,paste0("de-novosiganturesCellLinesTest",num.sign,".png")),width = 1000,height = 500)
# msigPlot(cell.signature.exp, main = "cell.lines")+
#   scale_fill_manual(values = c(col_vector[1:30]))
# 
# dev.off()
# 
# 
# #Another way to represent the cell.signature.exp results
# cse <- mutSignatures::as.data.frame(cell.signature.exp)
# df.sig <- c()
# for(i in seq_len(ncol(cse))){
#   sn <- colnames(cse)[i]
#   df.sig<- c(df.sig,cse[,sn])
#   
# }
# df.sig <- data.frame(Freq.Sig = df.sig, Signatures= rep(rownames(cse),8))
# df.sig$Cell.lines <- c(rep("S462",num.sign),
#                        rep("ST88-14",num.sign),
#                        rep("90-8-TL",num.sign),
#                        rep("SNF96-2",num.sign),
#                        rep("NMS-2",num.sign),
#                        rep("STS-26T",num.sign),
#                        rep("HS-Sch-2",num.sign),
#                        rep("HS-PSS",num.sign))
# df.sig$Cell.lines <- factor(df.sig$Cell.lines, levels = unique(df.sig$Cell.lines))               
# 
# svg(filename = file.path(image.dir, paste0("de-novosiganturesCellLinesTest_",num.sign,".svg")), width = 13, height = 9)
# 
# ggplot(df.sig, aes(Cell.lines, Freq.Sig, fill = Signatures)) + 
#   geom_bar(stat="identity", position = "stack") + 
#   scale_fill_brewer(palette = "Set2")+
#   theme(axis.text.x = element_text(angle = 90, size = 20, hjust = 1, vjust = 0),
#         axis.title = element_text(size = 30),
#         axis.text.y = element_text(size = 20),
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 30),
#         panel.background = element_rect(fill = "white", color="white", size = 2),
#         panel.grid= element_line(colour = "grey", size = 0.5, linetype = "solid") )
# dev.off()
# 
# 
# 
# msig1 <- matchSignatures(mutSign = cell.signature.sig,
#                          reference = cosmic, 
#                          threshold = 0.45,
#                          plot = T)
# png(file.path(image.dir,paste0("de-novoSigantures_MatchCosmic",num.sign, ".png")),width = 1000,height = 500)
# msig1
# dev.off()


# After obtaining the match signatures results,
# # we select those cosmic signatures with sum coefficient Beta >0.2.
# cell.sigs.cosmic <- mutSignatures::resolveMutSignatures(mutCountData = cells.signature.counts,
#                                                         signFreqData=cosmic)
# coef.beta <- mutSignatures::as.data.frame(cell.sigs.cosmic$results$freq.result)
# cosmic.selected<- rownames(coef.beta)[rowSums(coef.beta)>0.2]
# coef.beta <- coef.beta[cosmic.selected,]
# 
# png(file.path(image.dir,"Pheatmap_0.2CosmicSigantures_CoeffBeta.png"),width = 1000,height = 500)
# pheatmap(coef.beta)
# dev.off()
# # 
# cosmic.selected <- c("COSMIC.2","COSMIC.5","COSMIC.7","COSMIC.30")
# # cosmic.selected <- rownames(msig1$distanceMatrix)[rowSums(msig1$distanceMatrix)<0.6]
# cosmic.df <- mutSignatures::as.data.frame(cosmic)
# cosmic.df <- cosmic.df[,cosmic.selected]
# cosmic.selected <- mutSignatures::as.mutation.signatures(cosmic.df)
# # cosmic.selected <- cosmic

cell.sigs.cosmic <- mutSignatures::resolveMutSignatures(mutCountData = cells.signature.counts,
                                                        signFreqData=cosmic)

counts.cell.sig <- cell.sigs.cosmic$results$count.result

png(file.path(image.dir,"SelectedCosmicsigantures_CellLines.png"),width = 1000,height = 500)
msigPlot(counts.cell.sig, sample =colnames(counts.cell.sig) )+ 
  scale_fill_manual(values = c(col_vector[1:30]))
dev.off()


# coef.beta <- mutSignatures::as.data.frame(cell.sigs.cosmic$results$freq.result)
# pheatmap(coef.beta)
# 
# signatures <- as.numeric(gsub("COSMIC.","", colnames(as.data.frame(cosmic.selected))))
# for(i in seq_len(length(signatures))){
#  
#   png(file.path(image.dir,paste0("COSMIC_", signatures[i],"_SelectedCosmicsigantures.png")),width = 1000,height = 500)
#   
#   msigPlot(cosmic.selected, signature = i, ylim = c(0, 0.3))
#   dev.off()
# }
# 
# 
# msig2 <- matchSignatures(mutSign = cell.signature.sig,
#                          reference = cosmic.selected, 
#                          threshold = 0.45,
#                          plot = T)
# png(file.path(image.dir,"de-novoSigantures_MatchSelectedCosmic.png"),width = 1000,height = 500)
# msig2
# dev.off()



cse <- mutSignatures::as.data.frame(counts.cell.sig)
#Filtering those signatures without counts
cse <- cse[-which(rowSums(cse)==0),]

df.sig <- c()
for(i in seq_len(ncol(cse))){
  sn <- colnames(cse)[i]
  df.sig<- c(df.sig,cse[,sn])
  
}
df.sig <- data.frame("Mutation counts by signature" = df.sig, COSMIC.Signatures= rep(rownames(cse),8))
df.sig$COSMIC.Signatures <- gsub("COSMIC.","",df.sig$COSMIC.Signatures)
l.cosmic <- length(unique(df.sig$COSMIC.Signatures))
df.sig$Cell.lines <- c(rep("S462",l.cosmic),
                       rep("ST88-14",l.cosmic),
                       rep("90-8-TL",l.cosmic),
                       rep("SNF96-2",l.cosmic),
                       rep("NMS-2",l.cosmic),
                       rep("STS-26T",l.cosmic),
                       rep("HS-Sch-2",l.cosmic),
                       rep("HS-PSS",l.cosmic))

df.sig$Cell.lines <- factor(df.sig$Cell.lines, levels = unique(df.sig$Cell.lines)) 
df.sig$COSMIC.Signatures <- factor(df.sig$COSMIC.Signatures,levels= unique(df.sig$COSMIC.Signatures))

colnames(df.sig)
# df.sig$Mutation.counts.by.signature
set.seed(1214)
test.col <- sample(col_vector,nrow(cse))
test.col[1] <- "#6A9972"
test.col[6] <-"#8A0149"
test.col[15] <- "#595959"
test.col[17] <-"#943CF3"
# svg(filename = file.path(image.dir, paste0("de-novosigantures_MatchSelectedCosmic_",num.sign,".svg")), width = 13, height = 9)
svg(filename = file.path(image.dir, paste0("Fitting_MatchAllCosmic_",num.sign,".3.svg")), width = 13, height = 9)
ggplot(df.sig, aes(Cell.lines, Mutation.counts.by.signature, fill = COSMIC.Signatures)) +
  geom_bar(stat="identity", position = "stack") + 
  # scale_fill_brewer(palette = "primary.colors")+
  scale_fill_manual(values = test.col)+
  theme(axis.text.x = element_text(angle = 90, size = 20, hjust = 1, vjust = 0, color = "black"),
        axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 20,colour = "black"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 30),
        panel.background = element_rect(fill = "white", color="white", size = 2),
        panel.grid= element_line(colour = "grey", size = 0.5, linetype = "solid") )
dev.off()

