##################### Mut Signatures #############
library(mutSignatures)
# library(ggplot2)
# library(RColorBrewer)
# library(pheatmap)

execution.dir<- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_cellLines/WGS"
load(file.path(execution.dir,"results", "MPNSTcells_MutSignatureCounts.RData"))

# df<- mutSignatures::as.data.frame(cells.signature.counts)
# how many signatures should we extract? 
# num.sign <- "{%NUM_SIGN%}"
n.signatures <- c(7,5,2)
i=1
for(i in seq_len(length(n.signatures))){
  # Define parameters for the non-negative matrix factorization procedure.
  # you should parallelize if possible
  num.sign <- n.signatures[i]
  mutsig.dir <- file.path(execution.dir,"results","MutSignatures")
  if(!file.exists(mutsig.dir))dir.create(mutsig.dir)
  params <- 
    mutSignatures::setMutClusterParams( 
      num_processesToExtract = num.sign,    # num signatures to extract
      num_totIterations = 1000,               # bootstrapping: usually 500-1000
      num_parallelCores = 20)                # total num of cores to use (parallelization)
  
  # Extract new signatures - may take a while
  cell.signature.analysis <- 
    decipherMutationalProcesses(input = cells.signature.counts,
                                params = params)
  
  save(cell.signature.analysis, file = file.path(mutsig.dir,
                                                 paste0("MPNSTcells_MutSignatureAnalysis_1000Iter_", num.sign,"sign.RData")))
  
}
load(file.path(execution.dir,"results","MutSignatures",
               paste0("MPNSTcells_MutSignatureAnalysis_1000Iter_",num.sign,"sign.RData")))
