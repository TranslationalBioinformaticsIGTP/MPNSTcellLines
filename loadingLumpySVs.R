#MPNSTcellLines
# vcf.file <-file.path(execution.dir,"results/Lumpy_smoove_hg38/MPNST_cellLines-smoove.genotyped.vcf.gz")
# genome <- "hg38"
# sv.type <- "all"
# sample.name <- "S462"

loadingLumpySVs <- function(vcf.file, genome, sv.type, sample.name, verbose = TRUE){
  
  #test if vcf.file is a character
  if(is.character(vcf.file)){
    if(file.exists(vcf.file)){
      vcf.data <- VariantAnnotation::readVcf(vcf.file, genome = genome)
      
    }else{
      stop("file does not exist")
    }
    
  }else{
    stop("vcf.file parameter is not a character")
  }
  
  
  #test if sv.type is a character
  if(is.character(sv.type)){
    if(!sv.type %in% c("DUP","DEL","INV","ICR", "all")){
      stop("sv.type parameter must be either DUP, DEL, INV, ICR or all")
    }
  }else{
    stop("sv.type should be a character")
  }
  
  #Loading info and geno information from vcf
  info.lumpy <- info(vcf.data)
  geno.lumpy <- geno (vcf.data)
  geno.lumpy.GT <- geno.lumpy$GT
  
  
  rowranges.vcf <- rowRanges(vcf.data, genome = genome)
  rowranges.vcf <- UCSCStyle(rowranges.vcf)
  
  mcols(rowranges.vcf) <- mcols(rowranges.vcf)[,c(2:3)]
  mcols(rowranges.vcf) <- cbind(mcols(rowranges.vcf),
                                  IMPRECISE = info.lumpy$IMPRECISE,
                                  SVTYPE = info.lumpy$SVTYPE,
                                  SVLEN = info.lumpy$SVLEN, 
                                  END = info.lumpy$END, 
                                  MATEID = info.lumpy$MATEID,
                                  STRANDS = info.lumpy$STRANDS)
  
  #Filtering chromosomes
  rowranges.vcf <- regioneR::filterChromosomes(rowranges.vcf, organism = "hg")
  
  #Genotype = 0/0 means that there is no SV in that sample.
  #If SV is 0/1 means that is heterozygous 1/1 means homozygous SV.
  
  geno.lumpy.GT.s1 <- geno.lumpy.GT[,sample.name][!grepl(pattern = "0/0", x = geno.lumpy.GT[,sample.name])]
  geno.lumpy.GT.s1 <- geno.lumpy.GT.s1[!grepl(pattern = "[.]", x = geno.lumpy.GT.s1)]
  
  gt.names <- names(geno.lumpy.GT.s1)
  lumpy.sample <- rowranges.vcf[names(rowranges.vcf) %in% gt.names]
  
  if(verbose) message("Getting ", sv.type, " SVs...")
  lumpy.svs <- GRanges()
  
  # Duplications
  if(sv.type == "DUP"){

    lumpy.svs <- lumpy.sample[grepl(pattern = sv.type, x = lumpy.sample$SVTYPE)]
    
    #filtering out imprecise SVs
    # lumpy.svs <- lumpy.svs[lumpy.svs$IMPRECISE==FALSE]
    end(lumpy.svs) <- lumpy.svs$END
    if(length(lumpy.svs)>0) lumpy.svs$Sample <- sample.name

  }
  
  # Deletions
  if(sv.type == "DEL"){

    lumpy.svs <- lumpy.sample[grepl(pattern = sv.type, x = lumpy.sample$SVTYPE)]
    
    #filtering out imprecise SVs
    # lumpy.svs <- lumpy.svs[lumpy.svs$IMPRECISE==FALSE]
    end(lumpy.svs) <- lumpy.svs$END
    if(length(lumpy.svs)>0) lumpy.svs$Sample <- sample.name

  }

  # Translocations or Inv
  transloc.inv <- lumpy.sample[grepl(pattern = "BND", x = lumpy.sample$SVTYPE)]
  chromosomes <- seqlevels(lumpy.sample)
    
    for(chr in chromosomes){
      #chromosome by chromosome
      
      trans.inv.chr <- transloc.inv[as.character(seqnames(transloc.inv)) == chr]
      if(length(trans.inv.chr)== 0)next
     
       if(sv.type == "ICR"){
        # We select those chr different from the selected one.
        trans.chr <- trans.inv.chr[!grepl(pattern = paste0(chr,":"), x=unlist(trans.inv.chr$ALT))]
        if(length(trans.chr)== 0)next
        
        #filtering out imprecise SVs
        # trans.chr <- trans.chr[trans.chr$IMPRECISE =="FALSE"]
        
        #We proceed to find the chromosome where the translocation is produced
        transloc <- lumpy.sample[names(lumpy.sample) %in% as.character(unlist(trans.chr$MATEID))]
        transloc.position<- paste0(seqnames(transloc),":", start(transloc),"-",start(transloc))
        names(transloc.position) <- names(transloc)
        transloc.position <- transloc.position[as.character(unlist(trans.chr$MATEID))]
        trans.chr$END <- transloc.position
        lumpy.svs <- c(lumpy.svs,trans.chr)
        if(length(lumpy.svs)>0) lumpy.svs$Sample <- sample.name
      } 
      
      if(sv.type == "INV"){
        # We select the same chromosomes to detect inversions 
        inv <- trans.inv.chr[grepl(pattern = paste0(chr,":"), x= unlist(trans.inv.chr$ALT))]
        
        if(length(inv)== 0)next
        
        #filtering out imprecise SVs
        # inv <- inv[inv$IMPRECISE =="FALSE"]
        inv.end <- unlist(inv$ALT)
        
        inv.end <- lumpy.sample[names(lumpy.sample) %in% as.character(unlist(inv$MATEID))]
        inv.position<- paste0(seqnames(inv.end),":", start(inv.end),"-",start(inv.end))
        names(inv.position) <- names(inv.end)
        inv.position <- inv.position[as.character(unlist(inv$MATEID))]
        inv$END <- inv.position
        lumpy.svs <- c(lumpy.svs, inv)
        if(length(lumpy.svs)>0) lumpy.svs$Sample <- sample.name
      }
      
    }
   
    if(sv.type == "all"){
      
      lumpy.svs <- list()
      lumpy.svs[["DUP"]] <- loadingLumpySVs(vcf.file = vcf.file, genome = genome, sv.type = "DUP", sample.name = sample.name, verbose = FALSE)
      lumpy.svs[["DEL"]] <- loadingLumpySVs(vcf.file = vcf.file, genome = genome, sv.type = "DEL", sample.name = sample.name, verbose = FALSE)
      lumpy.svs[["ICR"]] <- loadingLumpySVs(vcf.file = vcf.file, genome = genome, sv.type = "ICR", sample.name = sample.name, verbose = FALSE)
      lumpy.svs[["INV"]] <- loadingLumpySVs(vcf.file = vcf.file, genome = genome, sv.type = "INV", sample.name = sample.name, verbose = FALSE)
      
      
      
    }
    message(sv.type, " SVs obtained")
    
  return(lumpy.svs)
  
}
