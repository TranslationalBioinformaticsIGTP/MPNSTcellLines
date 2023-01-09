#########################################
#       Preparation of Circos files     #
#########################################

######################## Functions ###########################

#If an error occurs, clean up as much as possible and exit with an error

changeTemplate <- function(x, template, replacement) {
  #message(paste0("Changing template   ", template, "   with   ", replacement ))
  x <- gsub(x=x, pattern=template, replacement=replacement, fixed=TRUE)
  return(x)
}
# First we need to generate a file per sample containing all ICR and INV
#INV and ICR will have an specific color.
execution.dir <- "./MPNST_cellLines/WGS"

#Creating circos directory
for(i in seq_len(length(sample.names))){
  circos.dir <- file.path(execution.dir,"results",sample.names[i],"Circos")
  if(!file.exists(circos.dir))dir.create(circos.dir)
}

load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.v4.RData")))
load(file= file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.v3.RData"))

# load(file = file.path(execution.dir,"results/",paste0("All_Lumpy_SV_Regions_FilteredCommonSVinSamples.RData")))
# load(file= file.path(execution.dir,"results/SVsCliffAndLumpy_AllSamples.RData"))

# filt.list$S462$AllRegions$SVTYPE
# filt.list$S462$AllRegionsFiltered$SVTYPE
# filt.list$S462$AllRegionsFiltered[filt.list$S462$AllRegionsFiltered$SVTYPE %in% c("ICR","INV")]
# join.cliff.lumpy$

circos.frst.data <-list()  
i=1
icr.col <- "color=(83,70,219)"
inv.col <- "color=(219,142,69)"

for(i in seq_len(length(filt.list))){
  sn <- names(filt.list)[i]
  lumpy.icr.inv <- filt.list[[sn]]$AllRegionsFiltered
  lumpy.icr.inv <- lumpy.icr.inv[lumpy.icr.inv$SVTYPE %in%c("ICR","INV")]
  mcols(lumpy.icr.inv) <-mcols(lumpy.icr.inv)[c(6,4,13)]
  cliff <- join.cliff.lumpy[[sn]]
  cliff <- cliff[cliff$SVTYPE%in%c("ICR","INV")]
  mcols(cliff) <- mcols(cliff)[c(1,2,3)]
  circos.data <- sort(unique(c(lumpy.icr.inv,cliff)))
  
  circos.data$SVTYPE[circos.data$SVTYPE == "ICR"] <- icr.col
  circos.data$SVTYPE[circos.data$SVTYPE == "INV"] <- inv.col
  circos.data.start <- data.frame(chr=as.character(seqnames(circos.data)), start= start(circos.data),end =end(circos.data))
  circos.data.start$chr <- gsub("chr","hs",circos.data.start$chr)
  circos.data.end <- toDataframe(toGRanges(circos.data$END))
  circos.data.end[,1] <- gsub("chr","hs",circos.data.end[,1])
  circos.data <- cbind(circos.data.start, circos.data.end,circos.data$SVTYPE)
  write.table(circos.data,file = file.path(execution.dir,"Circos_files",paste0("translocations_",sn,".txt")), sep="\t", col.names = FALSE,row.names = FALSE, quote = FALSE)
  
  
}

########## Modification of circos.transloc.conf #######
test <- "circos.transloc"
for (i in seq_len(length(sample.names))){
  sample.name <-sample.names[i]
  # cluster.file <- file.path(execution.dir, "toRunInCluster", paste0("Run_in_cluster_",test, ".sh"))
  cluster.file <- file.path(execution.dir,"Circos_files/WGS_MPNST_CL", paste0(test, ".conf"))
  cluster.file <-'./MPNST_cellLines/WGS/Circos_files/WGS_MPNST_CL/circos.transloc.conf'
  # t <- thr[sample.name,]
  # fq.dir <- sample.data$file.fastq[i]
  content <- readLines(cluster.file)
  # content <- changeTemplate(content, "{%SAMPLE_DIR}", execution.dir)
  content <- changeTemplate(content,"{%SAMPLE_NAME%}", sample.name)
  new.cluster.file <- file.path("./MPNST_cellLines/WGS/Circos_files/WGS_MPNST_CL", paste0(test,"_", sample.name,".conf"))
  
  write(content, file = new.cluster.file)
  
  ### Lunch pipeline
  # if (test == "StrelkaGerm"){
  #   # if (length( list.files(file.path("./MPNST_Tumors/WGS/results/", sample.name,"/Strelka_GermVariants"))) != 0) next
  # }else{
  #   if(length(list.files(file.path("/imppc/labs/eslab/mmagallon/Projects/Locus_CDKN2A/MPNST_Tumors/WGS/results/", sample.name,"/Manta_ManConfig_hg38"))) != 0) next
  # }
  
  # system(command = paste0("runcircos --conf=",new.cluster.file),wait = TRUE)
  
}