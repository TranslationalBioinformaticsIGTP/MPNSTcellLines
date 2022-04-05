#             Integrative Biology of MPNST Functions                #
#####################################################################
# #PArameters for testing the functions
# info.param <- c("Gene.refGene",
#                 "ExonicFunc.refGene",
#                 "Func.refGene",
#                 "avsnp150",
#                 "AAChange.refGene",
#                 "CLNALLELEID",
#                 "CLNDN",
#                 "CLNDISDB",
#                 "CLNREVSTAT",
#                 "CLNSIG",
#                 "SIFT_pred",
#                 "Polyphen2_HDIV_pred",
#                 "Polyphen2_HVAR_pred",
#                 "LRT_pred",
#                 "FATHMM_pred",
#                 "MutationTaster_pred",
#                 "MutationAssessor_pred",
#                 "InterVar_automated",
#                 "AF",
#                 "AF_popmax","GeneDetail.refGene")
# 
# 
# geno.param <- c("AD","DP")

############### getVAFStrelka2Germline ######################
# getVAF from Strelka2 germline calling
getVAFStrelka2Germline <- function(vcf.file){
  
  if(class(vcf.file) == "CollapsedVCF"){
    geno.info <- geno(vcf.file)
    
  }else{
    stop("vcf.file parameter should be a CollapsedVCF object")
  }
  
  # VAF calculation: VAF = AD / DP = Depth of Alt Allele / Total Depth 
  if(!"AD" %in% names(geno.info)){
    stop("AD is not in geno from vcf.file")
  }else{
    # DP <- geno.info$DP # total depth
    # cn <- colnames(DP)
    # DP <- geno.info$DP[,cn]
    # AD <- geno.info$AD #Depth of Alt Allele
    # cn <- colnames(AD)
    # AD <- AD[,cn]
    # AD <- sapply(AD,function(x){x[2]})
    # 
    # # VAF Computation 
    # vaf <-AD/DP
    ad <- geno.info$AD
    ad.ref <- unlist(lapply(ad, "[", 1))
    ad.alt <- unlist(lapply(ad, "[", 2))
    vaf <- ad.alt/(ad.alt + ad.ref)
    vaf[is.nan(vaf)] <- NA
    return(vaf)
    
  }
  
}

############### getVAFStrelka2Somatic ######################

getVAFStrelka2Somatic <- function(gr.vcf, variant.type, verbose = TRUE){
  #Adding VAF to the sample info
  #Gettin VAF annotation
  
  # if(class(vcf.file) == "CollapsedVCF"){
  #   geno.tumor <- geno(vcf.file)
  # 
  #   stop("vcf.file parameter should be a CollapsedVCF object")
  # }
  
  if(class(gr.vcf)=="GRanges"){
    if(variant.type =="snv"){
      strelka.som.vaf <- data.frame(mcols(gr.vcf)[colnames(mcols(gr.vcf))%in% c("AU.T","AU.N","TU.T","TU.N","GU.T","GU.N","CU.T","CU.N")])
      if(ncol(strelka.som.vaf)!=8){
        stop("gr.vcf containing somatic snvs must contain the following columns: AU.T,AU.N,TU.T,TU.N,GU.T,GU.N,CU.T,CU.N")
      }
      
    }else if(variant.type == "indels"){
      strelka.som.vaf <- data.frame(mcols(gr.vcf)[colnames(mcols(gr.vcf))%in% c("TAR.T","TAR.N","TIR.T","TIR.N")])
      if(ncol(strelka.som.vaf)!=4){
        stop("gr.vcf containing somatic indels must contain the following columns: TAR.T,TAR.N,TIR.T,TIR.N")
      }
    }
  }else{
    stop("gr.vcf must be a GRanges object")
    
  }
  
  if(is.character(variant.type)){
    variant.type <- variant.type
  }else{
    stop("variant.type must be a character, either snv, indels")
  }
  # if(!is.null(sel.pos)){
  #   if(!is.character(sel.pos)){
  #     stop("sel.pos must  be NULL or a character vector containing rownames of the variants")
  #   }
  # }
  
  
  
  
  #stracting ref and alt nucleotides
  ref.allele <- rownames(strelka.som.vaf)
  ref.allele <- strsplit(x= ref.allele,"_")
  ref.allele <- lapply(ref.allele,function(x) x[2])
  ref.allele <- lapply(ref.allele,function(x) strsplit(x,"/"))
  alt.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [2])))
  ref.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [1])))
  strelka.som.vaf$REF <- ref.allele
  strelka.som.vaf$ALT <- alt.allele
  
  if(verbose) message("Extracting VAF from somatic ", variant.type)
  
  
  if(variant.type == "snv"){
    
    for(i in seq_len(nrow(strelka.som.vaf))){
      n.vaf <- rownames(strelka.som.vaf)[i]
      ref <- strelka.som.vaf$REF[i]
      alt <- strelka.som.vaf$ALT[i]
      
      message("row ", i, " out of ", nrow(strelka.som.vaf))
      #VAF normal sample
      norm.counts <- strelka.som.vaf[i,grepl(".N", colnames(strelka.som.vaf))]
      ref.norm.counts <- norm.counts[grepl(ref,colnames(norm.counts))][,1]
      
      if(alt == "T"){
        alt.norm.counts <- norm.counts[grepl(paste0(alt,"U"),colnames(norm.counts))]
        
      }else{
        alt.norm.counts <- norm.counts[grepl(alt, colnames(norm.counts))][,1]
        
      }
      
      strelka.som.vaf$norm.VAF[i] <- as.numeric(alt.norm.counts/(ref.norm.counts+alt.norm.counts))
      
      
      #VAF tumor sample (somatic VAF)
      tumor.counts <- strelka.som.vaf[n.vaf,grepl("\\.T", colnames(strelka.som.vaf))]
      
      if(alt == "T" || ref == "T"){
        alt.tumor.counts <- tumor.counts[grepl(paste0(alt,"U"), colnames(tumor.counts))][,1]
        ref.tumor.counts <- tumor.counts[grepl(paste0(ref,"U"), colnames(tumor.counts))][,1]
        
      }else{
        alt.tumor.counts <- tumor.counts[grepl(alt, colnames(tumor.counts))][,1]
        ref.tumor.counts <- tumor.counts[grepl(ref,colnames(tumor.counts))][,1]
      }
      
      
      strelka.som.vaf$som.VAF[i] <- as.numeric(alt.tumor.counts/(ref.tumor.counts+alt.tumor.counts))
      
    }
    
  }else if(variant.type== "indels"){
    
    #VAF    
    strelka.som.vaf$norm.VAF <- strelka.som.vaf$TIR.N/(strelka.som.vaf$TIR.N+strelka.som.vaf$TAR.N)
    strelka.som.vaf$som.VAF <- strelka.som.vaf$TIR.T/(strelka.som.vaf$TIR.T+strelka.som.vaf$TAR.T)
    
  }
  
  
  return(strelka.som.vaf)
  
}

# test <- getVAFStrelka2Somatic(vcf.file = vcf, variant.type= "snv")
# test$som.VAF

# vcf.file.path <- strelka.res.path

############## loadingStrelka2AnnotVcf ######################
#Loading Strelka2 vcf functions for mutsignatures or for obtaining pathogenic variants
loadingStrelka2AnnotVcf <- function(vcf.file.path, wes.wgs, genome, geno.param, info.param, analysis.type, add.vaf=FALSE, verbose =TRUE){
  
  if(is.character(vcf.file.path)){
    if(!file.exists(vcf.file.path)){
      stop("vcf.path does not exist")
    }
  }else{
    stop("vcf.path must be a character")
  }
  
  if(!is.character(genome)){
    stop("genome parameter must be a character")
  }
  #TODO check gero.param info.param....
  # if(is.null(geno.param)){
  #   geno.param <- character()
  # }
  # 
  # if(is.null(info.param)){
  #   info.param <- character()
  # }
  
  
  if(wes.wgs == "wes"){
    if(analysis.type == "pathogenic"){
      # We select only those variants present in exons or affecting splicing and we delete the synonymous_SNVs
      filt <- FilterRules(list(Func.refGene = function(x) as.character(unlist(info(x)$Func.refGene)) %in% c("splicing","exonic"),
                               ExonicFunc.refGene = function(x) as.character(unlist(info(x)$ExonicFunc.refGene))%in% c(".","frameshift_insertion","frameshift_deletion", "stopgain", "starloss","nonsynonymous_SNV","nonframeshift_deletion","nonframeshift_insertion")))
      
    }else if (analysis.type =="mutsignatures"){
      filt <- FilterRules(list(Func.refGene = function(x) as.character(unlist(info(x)$Func.refGene)) %in% c("splicing","exonic")))
      
    }
    
    #Preparing the filter (file to load)
    filt2 <- filterVcf(vcf.file.path, genome = genome, tempfile() , param =ScanVcfParam(fixed = NA, geno = geno.param, info=info.param) , filters= filt, verbose = FALSE)
    
    
  }else if (wes.wgs == "wgs"){
    if(analysis.type == "pathogenic"){
      filt <- FilterRules(list(
        Func.refGene = function(x) as.character(unlist(info(x)$Func.refGene)) %in% c("splicing","exonic"), 
        ExonicFunc.refGene = function(x) as.character(unlist(info(x)$ExonicFunc.refGene))%in% c(".","frameshift_insertion","frameshift_deletion", "stopgain", "starloss","nonsynonymous_SNV","nonframeshift_deletion","nonframeshift_insertion")))
      
      filt2 <- filterVcf(vcf.file.path, genome = genome, tempfile() , param =ScanVcfParam(fixed = NA, geno = geno.param, info=info.param) , filters= filt, verbose = FALSE)
      
    }else if(analysis.type == "mutsignatures"){
      filt2   <-  vcf.file.path    
    }
    
    
  }
  #Loading strelka variants
  if(verbose)message("Loading strelka variants from ", wes.wgs, " for ", analysis.type, " analysis")
  vcf <- VariantAnnotation::readVcf(file = filt2, genome = genome, param = ScanVcfParam(fixed = NA, geno =geno.param, info=info.param), verbose = FALSE)
  
  #Remove temp file filt2 if wes data loaded
  if(wes.wgs == "wes"){
    file.remove(filt2)
  }
  
  
  #Transforming the data to GRanges
  gr.vcf <- rowRanges(vcf, genome = genome)
  mcols(gr.vcf) <- mcols(gr.vcf)[-1]
  gr.vcf$REF <- as.character(gr.vcf$REF)
  
  #Obtaining ALT columns
  ref.allele <- names(gr.vcf)
  ref.allele <- strsplit(x= ref.allele,"_")
  ref.allele <- lapply(ref.allele,function(x) x[2])
  ref.allele <- lapply(ref.allele,function(x) strsplit(x,"/"))
  alt.allele <- unlist(lapply(ref.allele,function(x) lapply(x,function(y)y [2])))
  gr.vcf$ALT <- as.character(alt.allele)
  
  
  #Adding info from vcf
  info.tumor <-  VariantAnnotation::info(vcf)
  mcols(gr.vcf) <- cbind(mcols(gr.vcf), info.tumor)
  
  #Adding VAF 
  if(add.vaf == TRUE){
    geno.info <- VariantAnnotation::geno(vcf)
    geno.info <- names(geno.info)
    if("AD" %in% geno.info){ 
      if(verbose)message("obtaining VAF")
      
      vaf <- getVAFStrelka2Germline(vcf)
      gr.vcf$VAF <- vaf
      
    }else{
      warning("we could not obtain VAF as AD was not in the vcf geno information")
    }
    
    mcols(gr.vcf) <- mcols(gr.vcf)[c(1,length(mcols(gr.vcf)),2:(length(mcols(gr.vcf))-1))]
    
  }
  
  
  #Unlisting mcols of granges
  for(j in seq_len(ncol(mcols(gr.vcf)))){
    
    if(class(mcols(gr.vcf)[[j]]) =="CompressedCharacterList"){
      ll <- lapply(mcols(gr.vcf)[[j]],function(x)length(x))
      mcols(gr.vcf)[[j]][ll >=2 ]<-lapply(mcols(gr.vcf)[[j]][ll>=2],function(x)paste0(x,collapse = ","))
      mcols(gr.vcf)[j] <- unlist(mcols(gr.vcf)[[j]])
      
    }
    
  }
  
  #Transforming AF_popmax value for filtering later
  if("AF_popmax"%in% colnames(mcols(gr.vcf))){
    suppressWarnings(gr.vcf$AF_popmax <-   as.numeric(gr.vcf$AF_popmax))
    gr.vcf$AF_popmax[is.na(gr.vcf$AF_popmax)] <- "."
  }
  
  
  # transforming CLNSIG column for make it readable
  if("CLNSIG" %in% colnames(mcols(gr.vcf))){
    ll <- lapply(gr.vcf$CLNSIG,function(x)length(x))
    gr.vcf$CLNSIG[ll >=2 ]<-lapply(gr.vcf$CLNSIG[ll>=2],function(x)paste0(x,collapse = ","))
    
  }
  
  message("Strelka2 annotated variants loaded from ", wes.wgs, " for ", analysis.type)
  gr.vcf <- filterChromosomes(gr.vcf)
  
  return(gr.vcf)
}

# vcf.file.path <-"/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WES/results/A09-121_cNF/Strelka_GermVariants2.9.10/results/variants/A09-121_cNF_annotated_PASS.vcf.gz"
# genome <- "hg38"
# wes.wgs <- "wes"
# analysis.type <- "mutsignatures"
# 
# 
# test <- loadingStrelka2AnnotVcf(vcf.file.path, wes.wgs, genome, geno.param = NULL, analysis.type = "mutsignatures",info.param = info.param )


############### loadingStrelka2SomaticVcf ######################
loadingStrelka2SomaticVcf <- function(vcf.file.path, variant.type, regions = NULL, add.vaf = FALSE, genome, verbose =TRUE){
  
  if(is.character(vcf.file.path)){
    if(!file.exists(vcf.file.path)){
      stop("vcf.path does not exist")
    }
  }else{
    stop("vcf.path must be a character")
  }
  
  if(!is.character(genome)){
    stop("genome parameter must be a character")
  }
  
  #TODO regions check
  if (!is.null(regions)) 
    regions <- tryCatch(regioneR::toGRanges(regions), error = function(e) {
      stop("regions must be in a valid format accepted by toGRanges.\n ", 
           e)
    })
  
  #Loading Strelka2 somatic variants
  if(is.null(regions)){
    vcf <- VariantAnnotation::readVcf(file = vcf.file.path, genome = genome, param = ScanVcfParam())
    
  }else{
    vcf <- VariantAnnotation::readVcf(file = vcf.file.path, genome = genome, param = ScanVcfParam(which = regions))
    
  }
  gr.vcf <- rowRanges(vcf, genome = genome)
  
  #getting VAF from somatic varinats
  #preparing data
  geno.tumor <- geno(vcf)
  if(variant.type == "snv"){
    strelka.som.vaf <- cbind(data.frame(AU.T = geno.tumor$AU[,,1][,"TUMOR"], AU.N = geno.tumor$AU[,,1][,"NORMAL"]),
                             data.frame(CU.T = geno.tumor$CU[,,1][,"TUMOR"],CU.N = geno.tumor$CU[,,1][,"NORMAL"]),
                             data.frame(GU.T = geno.tumor$GU[,,1][,"TUMOR"], GU.N = geno.tumor$GU[,,1][,"NORMAL"]),
                             data.frame(TU.T = geno.tumor$AU[,,1][,"TUMOR"], TU.N = geno.tumor$TU[,,1][,"NORMAL"]))
    
  }else if(variant.type == "indels"){
    if(nrow(data.frame(geno.tumor$TAR)<=1)){
      strelka.som.vaf <- cbind(data.frame(TAR.T = geno.tumor$TAR[,,1]["TUMOR"], TAR.N = geno.tumor$TAR[,,1]["NORMAL"])
                               ,data.frame(TIR.T = geno.tumor$TIR[,,1]["TUMOR"],TIR.N = geno.tumor$TIR[,,1]["NORMAL"]))
      
    }else{
      strelka.som.vaf <- cbind(data.frame(TAR.T = geno.tumor$TAR[,,1][,"TUMOR"], TAR.N = geno.tumor$TAR[,,1][,"NORMAL"])
                               ,data.frame(TIR.T = geno.tumor$TIR[,,1][,"TUMOR"],TIR.N = geno.tumor$TIR[,,1][,"NORMAL"]))
      
    }
    
  }
  
  
  mcols(gr.vcf) <- cbind(mcols(gr.vcf),strelka.som.vaf)
  
  if(add.vaf == TRUE){
    if(verbose)message("obtaining VAF from ", variant.type)
    vaf <-  getVAFStrelka2Somatic(gr.vcf = gr.vcf, variant.type = variant.type, verbose = verbose)
    mcols(gr.vcf) <- vaf[c("REF", "ALT", "norm.VAF", "som.VAF" )]
  }
  
  
  if(verbose)message("Strelka2 Somatic variants loaded from ", wes.wgs)
  
  return(gr.vcf)
}

# vcf.annotated.file <- vcf.file.path 


############### loadingAnnotStrelkaVCFforMutSig ######################

# Function for loading and proccesing vcf for mutSignatures from annovar annotated  calling
loadingAnnotStrelkaVCFforMutSig <- function(vcf.annotated.file, 
                                            genome,
                                            info.param,
                                            geno.param,
                                            wes.wgs,
                                            add.vaf=FALSE,
                                            verbose = TRUE){
  
  #check if vcf.annotated file is a character
  if(is.character(vcf.annotated.file)){
    if(file.exists(vcf.annotated.file)){
      gr.vcf <- loadingStrelka2AnnotVcf(vcf.file.path = vcf.annotated.file,
                                        wes.wgs = wes.wgs,
                                        genome = genome,
                                        geno.param = geno.param,
                                        info.param = info.param,
                                        add.vaf= add.vaf,
                                        analysis.type = "mutsignatures",
                                        verbose =verbose)
    }else{
      stop(vcf.annotated.file, "does not exist")
    }
    
  }
  
  if(verbose) message("Variants ready for filtering before running mutSignatures ")
  return(gr.vcf)
  
}

#####################loadingAnnotStrelkaVCFforPathogeneity#########
# Function for loading and proccesing vcf for mutSignatures from annovar annotated  calling
loadingAnnotStrelkaVCFforPathogeneity <- function(vcf.annotated.file, 
                                                  genome,
                                                  info.param,
                                                  geno.param,
                                                  wes.wgs,
                                                  add.vaf=FALSE,
                                                  verbose = TRUE){
  
  #check if vcf.annotated file is a character
  if(is.character(vcf.annotated.file)){
    if(file.exists(vcf.annotated.file)){
      gr.vcf <- loadingStrelka2AnnotVcf(vcf.file.path = vcf.annotated.file,
                                        wes.wgs = wes.wgs,
                                        genome = genome,
                                        geno.param = geno.param,
                                        info.param = info.param,
                                        analysis.type = "pathogenic",
                                        add.vaf=add.vaf,
                                        verbose =verbose)
    }else{
      stop(vcf.annotated.file, "does not exist")
    }
    
  }
  
  if(verbose) message("Splicing and exonic variants ready for filtering from ", wes.wgs, " data")
  return(gr.vcf)
  
}


############### loadingSomaticStrelkaVCFMutsig ######################
# vcf.annotated.file <- "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WES/results/GFS_MPNST/Strelka_SomaticVariants2.9.10/results/variants/GFS_MPNST_SNV_somaticVariants_annotated.vcf.gz"
# vcf.som.file <-  "/imppc/labs/eslab/mmagallon/Projects/Integrative_Biology/MPNST_tumors/WES/results/GFS_MPNST/Strelka_SomaticVariants2.9.10/results/variants/somatic.snvs_GFS_MPNST.vcf.gz"
# variant.type <- c("snv", "indels")
# variant.type <- "snv"

loadingSomaticStrelkaVCFforMutSig <- function(vcf.annotated.file,
                                              vcf.som.file, 
                                              genome,
                                              info.param,
                                              geno.param = NULL,
                                              variant.type,
                                              add.vaf = FALSE,
                                              wes.wgs,
                                              verbose = TRUE){
  
  #check if vcf.annotated file is a character
  if(is.character(vcf.annotated.file)){
    if(file.exists(vcf.annotated.file)){
      gr.vcf <- loadingStrelka2AnnotVcf(vcf.file.path = vcf.annotated.file,
                                        wes.wgs = wes.wgs,
                                        genome = genome,
                                        geno.param = geno.param,
                                        info.param = info.param,
                                        add.vaf = FALSE,
                                        analysis.type = "mutsignatures",
                                        verbose = verbose)
    }else{
      stop(vcf.annotated.file, "does not exist")
    }
    
  }
  
  if(is.character(vcf.som.file)){
    if(file.exists(vcf.som.file)){
      gr.vcf.som <- loadingStrelka2SomaticVcf(vcf.file.path = vcf.som.file,
                                              genome = genome,
                                              variant.type = variant.type,
                                              add.vaf= add.vaf,
                                              verbose = verbose)
    }else{
      stop(vcf.som.file, "does not exist")
    }
    
  }
  #Adding VAF to gr.vcf
  gr.vcf.som <- gr.vcf.som[names(gr.vcf.som)%in% names(gr.vcf)]
  mcols(gr.vcf) <- mcols(gr.vcf)[-1]
  mcols(gr.vcf) <- cbind(mcols(gr.vcf.som), mcols(gr.vcf))
  
  if(verbose) message("Variants ready for filtering before running mutSignatures ")
  
  return(gr.vcf)
  
}

# snv.som.mutsig <- loadingSomaticStrelkaVCFforMutSig(vcf.annotated.file = vcf.annotated.file,
#                                                     vcf.som.file = vcf.som.file,
#                                                     genome = genome,
#                                                     info.param = info.param,
#                                                     geno.param = NULL,
#                                                     variant.type = "snv",
#                                                     wes.wgs = "wes")

#####################loadingSomaticStrelkaVCFforPathogeneity#########
loadingSomStrelkaVCFforPathogeneity <- function(vcf.annotated.file,
                                                vcf.som.file, 
                                                genome,
                                                info.param,
                                                geno.param = NULL,
                                                variant.type,
                                                add.vaf = FALSE,
                                                wes.wgs,
                                                verbose = TRUE){
  
  #check if vcf.annotated file is a character
  if(is.character(vcf.annotated.file)){
    if(file.exists(vcf.annotated.file)){
      gr.vcf <- loadingStrelka2AnnotVcf(vcf.file.path = vcf.annotated.file,
                                        wes.wgs = wes.wgs,
                                        genome = genome,
                                        geno.param = geno.param,
                                        info.param = info.param,
                                        add.vaf = FALSE,
                                        analysis.type = "pathogenic",
                                        verbose =verbose)
    }else{
      stop(vcf.annotated.file, "does not exist")
    }
    
  }
  
  if(is.character(vcf.som.file)){
    if(file.exists(vcf.som.file)){
      gr.vcf.som <- loadingStrelka2SomaticVcf(vcf.file.path = vcf.som.file,
                                              genome = genome,
                                              variant.type = variant.type,
                                              regions = gr.vcf,
                                              add.vaf= add.vaf,
                                              verbose = verbose)
    }else{
      stop(vcf.som.file, "does not exist")
    }
    
  }
  #Adding VAF to gr.vcf
  gr.vcf.som <- gr.vcf.som[names(gr.vcf.som)%in% names(gr.vcf)]
  mcols(gr.vcf) <- mcols(gr.vcf)[-1]
  mcols(gr.vcf) <- cbind(mcols(gr.vcf.som), mcols(gr.vcf))
  
  if(verbose) message("Splicing and exonic variants ready for filtering from ", wes.wgs, " data")
  
  return(gr.vcf)
  
}
# snv.som.patho <- loadingSomaticStrelkaVCFforPathogeneity(vcf.annotated.file = vcf.annotated.file,
#                                                          vcf.som.file = vcf.som.file,
#                                                          genome = genome,
#                                                          info.param = info.param,
#                                                          geno.param = NULL,
#                                                          variant.type = "snv",
#                                                          wes.wgs = "wes")

###############computingPathogenicProp########################

# predictors <- info.param[10:18]
computingPredictorsPathogenicProp <- function(gr.vcf, predictors, verbose = TRUE){
  #Testing the parameters
  if(class(gr.vcf)!="GRanges"){
    stop("gr.vcf must be a GRanges object")
  }
  
  if(is.character(predictors)){
    if(FALSE %in%(predictors %in% colnames(mcols(gr.vcf)))){
      warning("not all given predictors are in the GRanges object")
    }
  }
  #Creating a data.frame for computing proportions
  df <- data.frame(mcols(gr.vcf)[predictors]) # selecting the introduced predictors
  df <- data.frame(apply(df,2,as.character))
  rownames(df) <- names(gr.vcf)
  
  if(verbose)("Obtaining pathogenic proportions")
  for(i in seq_len(ncol(df))){
    sub.df <- df[,i]
    patho.pos <- which(df[,i] %in%c("D","P","A","D","P","H","M"))
    
    if(length(patho.pos)==0)next
    
    df[patho.pos,i] <- "Pathogenic"
  }
  
  #Counting the number of pathogenic predictors
  if("CLNSIG" == colnames(df) || "InterVar_automated" == colnames(df)){
    df.patho <- df[,-which(colnames(df) %in%c("CLNSIG","InterVar_automated"))]# this are not predictors, they asign in other way
  }else{
    df.patho <- df
  }
  patho <- rowSums(df.patho == "Pathogenic") 
  uncert <- rowSums(df.patho == ".")
  prop <- (patho)/(ncol(df.patho)-uncert)
  
  gr.vcf$prop_pathogenic_pred <- prop
  
  if(verbose)("pathogenic proportion added to GRanges object")
  
  return(gr.vcf)
}

# ss <- computingPredictorsPathogenicProp(gr.vcf = snv.som.patho, predictors = predictors)






