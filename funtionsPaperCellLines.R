#Function kpPlotMeanCoveragebySample.
#This function was created for representing the mean of the number of alleles gained or lost
#We only add a parameter to kpPlotCoverage function from karyoploteR to easier represent the plot we wanted for
# cell lines characterization paper.

kpPlotMeanCoveragebySample <- function(karyoplot, data, num.samples, show.0.cov=TRUE, data.panel=1, r0=NULL, r1=NULL, col="#0e87eb", border=NULL, ymax=NULL, clipping=TRUE, ...) {
  #Check parameters
  #karyoplot
  if(missing(karyoplot)) stop("The parameter 'karyoplot' is required")
  if(!methods::is(karyoplot, "KaryoPlot")) stop("'karyoplot' must be a valid 'KaryoPlot' object")
  
  #data
  if(missing(data)) stop("The parameter 'data' is required")
  #TODO: If data is not a SimpleRleList, try to convert it to a GRanges before testing with is so other Region Set formats can be used
  
  if(!methods::is(data, "GRanges") && !methods::is(data, "SimpleRleList")) {
    data <- tryCatch(toGRanges(data), error=function(e) {stop("'data' must be a GRanges object or a SimpleRleList")})
  }  
  
  #Compute (if needed) the coverage
  #If its not a coverage object,  it's a GRanges. Compute the coverage
  if(!methods::is(data, "SimpleRleList")) { 
    #remove any region not in the currently used genome
    data <- data[seqnames(data) %in% karyoplot$chromosomes,]
    #Old version, problems when data had no seqinfo - data <- GenomeInfoDb::keepSeqlevels(data, karyoplot$chromosomes, pruning.mode="coarse")
    #Remove any unused seq level from the GRanges to fix problems with coverage and witdh
    seqlevels(data) <- karyoplot$chromosomes
    #the width parameter is needed so the coverage extends to the end of the chromosomes
    data <- GenomicRanges::coverage(data, width=karyoplot$chromosome.lengths[seqlevels(data)]) 
  }
  
  coverage.gr <- toGRanges(data)
  
  if(show.0.cov==FALSE) {
    coverage.gr <- coverage.gr[coverage.gr$coverage!=0]
  }
  
  if(is.null(ymax)) ymax <- max(max(coverage.gr$coverage))
  
  if(is.null(border)) border <- col
  
  # kpBars(karyoplot=karyoplot, data=coverage.gr,
  #        y0=0, y1=coverage.gr$coverage, ymin=0, ymax=ymax,
  #        r0=r0, r1=r1, data.panel=data.panel,
  #        col=col, border=border, clipping=clipping, ...)
  
  
  #To get kpArea to plot the real coverage (flat tops), we need to build a
  #GRanges with two elements per range, one at the start and one at the end
  #NOTE: this breaks the show.cov.0=FALSE. 
  #TODO: Make kpArea and kpLines deal with NAs in the data and set coverage=0 
  #to NA here
  cov.start <- coverage.gr
  end(cov.start) <- start(cov.start)
  cov.end <- coverage.gr
  start(cov.end) <- end(cov.end)
  cov.to.plot <- sort(c(cov.start, cov.end))
  
  #mean of coverage gain by sample
  cov.to.plot$coverage <- cov.to.plot$coverage/num.samples
  max(cov.to.plot$coverage)
  kpArea(karyoplot=karyoplot, data=cov.to.plot, y=cov.to.plot$coverage, 
         base.y = 0, ymin=0, ymax=ymax,
         r0=r0, r1=r1, data.panel=data.panel,
         col=col, border=border, clipping=clipping, ...)
  
  karyoplot$latest.plot <- list(funct="kpPlotCoverage", computed.values=list(max.coverage=max(max(coverage.gr$coverage)),
                                                                             coverage=coverage.gr,
                                                                             ymax=ymax))
  
  invisible(karyoplot)
}