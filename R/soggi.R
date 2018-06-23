#' The soggi function is the constructor for ChIPprofile objects.
#'
#' @name soggi
#' @rdname ChIPprofile
#' @export
#' @import methods reshape2 ggplot2 BiocGenerics S4Vectors IRanges GenomeInfoDb GenomicRanges Biostrings Rsamtools GenomicAlignments rtracklayer preprocessCore chipseq BiocParallel
#' @include allClasses.r plots.R peakTransforms.r
regionPlot <- function(bamFile,testRanges,samplename=NULL,nOfWindows=100,FragmentLength=150,style="point",distanceAround=NULL,distanceUp=NULL,distanceDown=NULL,distanceInRegionStart=NULL,distanceOutRegionStart=NULL,distanceInRegionEnd=NULL,distanceOutRegionEnd=NULL,paired=FALSE,normalize="RPM",plotBy="coverage",removeDup=FALSE,verbose=TRUE,format="bam",seqlengths=NULL,forceFragment=NULL,method="bin",genome=NULL,cutoff=80,downSample=NULL,minFragmentLength=NULL,maxFragmentLength=NULL){
  if(!verbose){
    suppressMessages(runRegionPlot())
  }
  result <- runRegionPlot(bamFile,testRanges,samplename,nOfWindows,FragmentLength,style,distanceAround,distanceUp,distanceDown,distanceInRegionStart,distanceOutRegionStart,distanceInRegionEnd,distanceOutRegionEnd,paired,normalize,plotBy,removeDup,format,seqlengths,forceFragment,method,genome,cutoff,downSample,minFragmentLength,maxFragmentLength)
  return(result)  
}

runRegionPlot <- function(bamFile,testRanges,samplename=NULL,nOfWindows=100,FragmentLength=150,style="point",distanceAround=NULL,distanceUp=NULL,distanceDown=NULL,distanceInRegionStart=NULL,distanceOutRegionStart=NULL,distanceInRegionEnd=NULL,distanceOutRegionEnd=NULL,paired=FALSE,normalize="RPM",plotBy="coverage",removeDup=FALSE,format="bam",seqlengths=NULL,forceFragment=NULL,method="bin",genome=NULL,cutoff=80,downSample=NULL,minFragmentLength=NULL,maxFragmentLength=NULL){

  # bamFile <- "~/Downloads/ENCFF049TYL.bam"
  # #bamFile <-"Downloads//mergedETOH.bwRange5.bw"
  # #bamFile <-"/Users/tcarroll//Downloads//Sample_R1-6hDupMarkedNormalised.bw"
  # library(GenomicRanges)
  # MelPeaks <- read.delim("~/Downloads/ENCFF591LSO.bed.gz",sep="\t",h=F)
  # MelGR <- GRanges(MelPeaks[,1],IRanges(MelPeaks[,2],MelPeaks[,3]))
  # library(Rsamtools)
  # # indexBam("~/Downloads/ENCFF049TYL.bam")
  # # indexBam("~/Downloads/ENCFF006JXP.bam")
  # # melSi <- regionPlot("~/Downloads/ENCFF049TYL.bam"
  # testRanges <- MelGR
  # nOfWindows=100
  # FragmentLength=150
  # style="point"
  # distanceAround=1500
  # distanceAround=20
  # distanceUp=NULL
  # distanceDown=NULL
  # distanceInRegionStart=1500
  # distanceOutRegionStart=1500
  # distanceInRegionEnd=1500
  # distanceOutRegionEnd=1500
  # paired=F
  # normalize="RPM"
  # plotBy="coverage"
  # removeDup=F
  # format="bam"
  # seqlengths=NULL
  if(format == "bam"){
    if(file.exists(bamFile) & is.na(index(BamFile(bamFile)))){
      message("Creating index for ",bamFile)
      indexBam(bamFile)
      message("..done")
    }
  }
  
  # testRanges <- soGGi:::GetGRanges(testRanges)
  ## Check parameters
  if(style != "percentOfRegion"){
    if(is.null(distanceAround)){
      distanceAround = 1500
    }
    if(is.null(distanceUp)){
      distanceUp <- distanceAround
    }
    if(is.null(distanceDown)){
      distanceDown <- distanceAround
    }
  }else{
    if(is.null(distanceAround)){
      distanceAround = 100
    }
    if(is.null(distanceUp)){
      distanceUp <- distanceAround
    }
    if(is.null(distanceDown)){
      distanceDown <- distanceAround
    }
  }
  
  if(is.null(distanceInRegionStart)){
    distanceInRegionStart=750
  }
  if(is.null(distanceOutRegionStart)){
    distanceOutRegionStart=1500
  }
  if(is.null(distanceInRegionEnd)){
    distanceInRegionStart=750
  }
  if(is.null(distanceOutRegionEnd)){
    distanceInRegionStart=1500
  }
  
  ## Initialize empty matrices and paramaters for collecting coverage analysis
  ## Find maximum distance to use for filtering out of bounds extended GRanges
  if(style == "region" | style=="regionandpoint"){
    posRegionStartMat <- NULL
    posRegionEndMat <- NULL
    negRegionStartMat <- NULL
    negRegionEndMat <- NULL
    RegionsMat <- NULL
    maxDistance <- max(distanceOutRegionStart,distanceOutRegionEnd)
    distanceUpStart <- distanceOutRegionStart
    distanceDownEnd <- distanceOutRegionEnd
    
  }
  
  if(style == "point"){
    PosRegionMat <- NULL
    NegRegionMat <- NULL
    RegionsMat <- NULL
    whatIsMax <- max(distanceAround,distanceUp,distanceDown)
    maxDistance <- whatIsMax
    distanceUpStart <- distanceUp
    distanceDownEnd <- distanceDown    
  }
  
  if(style == "percentOfRegion"){
    maxDistance <- round((distanceAround/100)*width(testRanges))
    RegionsMat <- NULL    
    distanceUpStart <- NULL
    distanceDownEnd <- NULL    
  }
  totalReads <- NA
  
  
  ## If format is bam, read header and get contig information
  
  if(format == "bam"){
    ## Get all chromosomes in bamFile
    message("Reading Bam header information...",appendLF = FALSE)
    allchrs <- names(scanBamHeader(bamFile)[[1]]$targets)
    lengths <- as.vector(scanBamHeader(bamFile)[[1]]$targets)
    names(lengths) <- allchrs
    message("..Done")
  }
  
  ## For remaining formats import data and find contig information
  
  # Import bigwig.
  
  # if(format=="bigwig"){
  #   message("Importing BigWig...",appendLF = FALSE)
  #   genomeCov <- import.bw(bamFile,as = "RleList")
  #   if(is.null(seqlengths)){
  #     seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
  #   }else{
  #     seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
  #   }
  #   lengths <- seqlengths(genomeCov)
  #   allchrs <- names(lengths)
  #   message("..Done")
  # }

  # bamFile <- "~/Downloads/ENCFF259VHD.bigWig"
  
  if(format=="bigwig"){
    message("Reading BigWig contig information...",appendLF = FALSE)
    bwFF <- BigWigFile(bamFile)
    # if(is.null(seqlengths)){
    #   seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
    # }else{
    #   seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
    # }

    lengths <- seqlengths(bwFF)
    allchrs <- names(lengths)
    message("..Done")
  }
  
  # If format is pwm PWM, calculate motifs on forward and reverse strand.
  
  if(format=="pwm"){
    bamFile <- pwmToCoverage(bamFile,genome,min=cutoff,removeRand=FALSE)
    format <- "rlelist"   
  }
  
  # If format is granges, simply covert to coverage (This would work for GenomicInterval or GenomicAlignments)
  
  if(format=="granges"){
    genomeCov <- coverage(bamFile)
    format <- "rlelist"   
  }  
  
  # If format is rlelist, import rle and set widths/contigs by seqlengths.
  
  if(format=="rlelist"){
    message("Importing rlelist",appendLF = FALSE)
    genomeCov <- bamFile
    if(is.null(seqlengths)){
      seqlengths(genomeCov) <- unlist(lapply(genomeCov,length))
    }else{
      seqlengths(genomeCov)[match(names(lengths),names(genomeCov))] <- lengths
    }
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..Done")
  }
  
  # Exclude and count regions which when extended are outside contig boundaries.
  
  if(style != "percentOfRegion"){
    ## Filter testRanges to those contained within chromosomes.
    message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
    testRangeNames <- unique(seqnames(testRanges))
    temptestranges <- GRanges()
    for(i in 1:length(testRangeNames)){
      perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
      temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+maxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                      & start(perchrRanges)-maxDistance > 0 ])
      #print(i)
    }
  }
  if(style == "percentOfRegion"){
    message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
    testRangeNames <- unique(seqnames(testRanges))
    temptestranges <- GRanges()
    for(i in 1:length(testRangeNames)){
      perChrMaxDistance <- maxDistance[as.vector(seqnames(testRanges) %in% as.vector(testRangeNames[i]))]
      perchrRanges <- testRanges[seqnames(testRanges) %in% as.vector(testRangeNames[i])]
      temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                      & start(perchrRanges)-perChrMaxDistance > 0 ])
      #print(i)
      perChrMaxDistance <- perChrMaxDistance[end(perchrRanges)+perChrMaxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                             & start(perchrRanges)-perChrMaxDistance > 0 ]
      distanceUpStart <- c(distanceUpStart,perChrMaxDistance)
    }
    distanceDownEnd <- distanceUpStart
    
  }  
  message("..Done")
  message("Filtered ",length(testRanges)-length(temptestranges)," of ",length(testRanges)," regions")
  testRanges <- temptestranges
  temptestranges <- NULL
  
  # Split ranges into +/- strand. Regions with no strand information are assigned to + strand
    
  message("Splitting regions by Watson and Crick strand..",appendLF = FALSE)
  mcols(testRanges) <- cbind(mcols(testRanges),data.frame(giID = paste0("giID",seq(1,length(testRanges)))))
  strand(testRanges[strand(testRanges) == "*"]) <- "+"
  testRangesPos <- testRanges[strand(testRanges) == "+"]
  testRangesNeg <- testRanges[strand(testRanges) == "-"]
  message("..Done")
  if(style=="percentOfRegion"){
    distanceUpStartPos <- distanceUpStart[as.vector(strand(testRanges) == "+")] 
    distanceDownEndPos <- distanceUpStartPos    
    distanceUpStartNeg <- distanceUpStart[as.vector(strand(testRanges) == "-")]
    distanceDownEndNeg <- distanceUpStartNeg
    message("..Done")
  }else{
    distanceUpStartPos <- distanceUpStart 
    distanceDownEndPos <- distanceDownEnd    
    distanceUpStartNeg <- distanceUpStart
    distanceDownEndNeg <- distanceDownEnd    
    message("..Done")
  }
  
#  if(style=="percentOfRegion"){
#    message("Filtering regions which are smaller than windows into region...",appendLF = FALSE)
#    ## Split Regions into those on positive and negative strands..
#    testRangesPos <- testRangesPos[width(testRangesPos) > nOfWindows]
#    testRangesNeg <- testRangesNeg[width(testRangesPos) > nOfWindows]
#    message("..Done")
#  }  
  
  if(style=="region"){
    message("Filtering regions which are smaller than windows into region...",appendLF = FALSE)
    ## Split Regions into those on positive and negative strands..
    testRangesPos <- testRangesPos[(end(testRangesPos)-distanceInRegionEnd) - (start(testRangesPos)+distanceInRegionStart) > nOfWindows]
    testRangesNeg <- testRangesNeg[(end(testRangesNeg)-distanceInRegionStart) - (start(testRangesNeg)+distanceInRegionEnd) > nOfWindows]
    message("..Done")
  }  
  
  message("Found ",length(testRangesPos)," Watson strand regions")
  message("Found ",length(testRangesNeg)," Crick strand regions")
  
  
  ## Extend regions and get positive versus negative regions
  message("Extending regions..",appendLF=FALSE)    
  exttestRanges <- c(GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStartPos,end(testRangesPos)+distanceDownEndPos)),
                     GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceDownEndNeg,end(testRangesNeg)+distanceUpStartNeg))
  )
  message("...done")   
  
  ## Create GRanges to be used in scanBamParam while reading in Bamfile regions.
  reducedExtTestRanges <- reduce(exttestRanges)

  ## Set up scanBanParam for reading in bam file.
  
  if(!removeDup){
    Param <- ScanBamParam(which=GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
  }else{
    Param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE),which=GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
  }
  
  ## if format is bam read in bamfile.
  if(format == "bam"){
    message("Reading tags from ",bamFile,appendLF=FALSE)
    #totalReads <- alignmentStats(bamFile)[,"mapped"]
    # QuasR is breaking builds so for now relies on USER input
    totalReads <- 10^6
  ##  if data is single end reads then import reads and reset reads to fragment length.
  ##  Calculate fragment length from cross-coverage if not provided
  
    if(paired==FALSE){
      total <- readGAlignments(bamFile,param=Param)
      message("..Done.\nRead in ",length(total)," reads")
      
      if(is.null(FragmentLength)){
        FragmentLength <- getShifts(total,lengths,shiftWindowStart=1,shiftWindowEnd=400)
      }
      
      message("Extending reads to fragmentlength of ",FragmentLength,appendLF=FALSE)
      temp <- resize(as(total,"GRanges"),FragmentLength,"start")
      message("..done")
    }

    ##  if data is paired end reads then import reads.
    ##  Reset fragment size and/or filter to fragment length range is specified.
  
    if(paired==TRUE){
      
      gaPaired <- readGAlignments(bamFile, 
                                         param=ScanBamParam(what=c("mpos"),
                                                            flag=scanBamFlag(isProperPair = TRUE,isFirstMateRead = TRUE),
                                                            mapqFilter=30))      
      tempPos <- GRanges(seqnames(gaPaired[strand(gaPaired) == "+"]),
                         IRanges(
                           start=start(gaPaired[strand(gaPaired) == "+"]),
                           end=mcols(gaPaired[strand(gaPaired) == "+"])$mpos
                           +qwidth(gaPaired[strand(gaPaired) == "+"])))
      tempNeg <- GRanges(seqnames(gaPaired[strand(gaPaired) == "-"]),
                         IRanges(
                           start=mcols(gaPaired[strand(gaPaired) == "-"])$mpos,                        
                           end=end(gaPaired[strand(gaPaired) == "-"])
                         )) 
      temp <- c(tempPos,tempNeg)
      message("..Done.\nRead in ",length(temp)," reads")
      if(removeDup){
        message("Removing duplicates")
        beforeDupR <- length(temp)
        temp <- unique(temp)
        AfterDupR <- length(temp)
        message("Removed ",beforeDupR-AfterDupR," duplicates")
      }
      #temp <- GRanges(seqnames(tempPaired),IRanges(start(left(tempPaired)),end(right(tempPaired))))

      if(!is.null(minFragmentLength)){
        temp <- temp[width(temp) > minFragmentLength]
      }
      if(!is.null(maxFragmentLength)){
        temp <- temp[width(temp) < maxFragmentLength]
      }   
      if(!is.null(forceFragment)){
        message("Forcing fragments to be centred and set to ",forceFragment,"..",appendLF=FALSE)
        temp <- resize(temp,forceFragment,"center")
        message("..done")        
      }

      message("..done")
    }
  
    ## Downsample single or paired end reads if specified and create coverage rlelist
    message("Calculating coverage..",appendLF=FALSE)
    
    seqlengths(temp)[match(names(lengths),names(seqlengths(temp)))] <- lengths
    if(!is.null(downSample)){
        if(downSample < 1 & downSample > 0){
          temp <- temp[sample(length(temp),round(length(temp))*downSample),]
        }else if(downSample > 1){
          temp <- temp[sample(length(temp),downSample),]
        } 
    }

    genomeCov <- coverage(temp)
    lengths <- seqlengths(genomeCov)
    allchrs <- names(lengths)
    message("..done")
  }
  
  if(format=="bigwig"){
    bwSelect <- BigWigSelection(GRanges(seqnames=seqnames(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),IRanges(start=start(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]),end=end(reducedExtTestRanges[seqnames(reducedExtTestRanges) %in% allchrs]))))
    genomeCov <- import.bw(bamFile,selection=bwSelect,as="RleList")
  }
  
  
  chromosomes <- seqlevels(genomeCov) 
  chromosomes <- chromosomes[chromosomes %in% unique(seqnames(reducedExtTestRanges))]
  
  # If style is "point" creates matrix of  per base pair coverage around centre of GRanges
  
  if(style=="point"){
    testRangesPos <- resize(testRangesPos,1,"center")
    testRangesNeg <- resize(testRangesNeg,1,"center")
    RangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceUpStart,start(testRangesPos)+distanceDownEnd),strand=Rle("+",length(testRangesPos)),mcols(testRangesPos))
    RangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(end(testRangesNeg)-distanceDownEnd,end(testRangesNeg)+distanceUpStart),strand=Rle("-",length(testRangesNeg)),mcols(testRangesNeg))  
    message("Calculating coverage across regions\nCalculating per contig. ")
    
    for(c in 1:length(chromosomes)){
      message(paste0("contig: ",c))      
      if(length(RangesPos[seqnames(RangesPos) %in% chromosomes[c]]) > 0){
        PosRegionMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(RangesPos[seqnames(RangesPos) %in% chromosomes[c]])]),ncol=mean(width(RangesPos)),byrow=TRUE)
        rownames(PosRegionMat) <- RangesPos[seqnames(RangesPos) %in% chromosomes[c]]$giID
      }
      if(length(RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]]) > 0){
        NegRegionMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]])])),ncol=mean(width(RangesNeg)),byrow=TRUE)
        rownames(NegRegionMat) <- RangesNeg[seqnames(RangesNeg) %in% chromosomes[c]]$giID
      }    
      RegionsMat <- rbind(RegionsMat,PosRegionMat,NegRegionMat)
      PosRegionMat <- NULL
      NegRegionMat <- NULL
    }
    message("Creating ChIPprofile.")
    
    ## Create matrix as summarisedexperiment with column and row names to be exported as part of ChIPprofile object
    require(DelayedMatrixStats)
    # profileMat <- realize(RegionsMat,"HDF5Array")
    profileMat <- RegionsMat
    colnames(profileMat) <- c(paste0("Point_Centre",seq(0-distanceUpStart,-1)),"Point_Centre",paste0("Point_Centre",seq(1,distanceDownEnd)))
    filteredRanges <- c(RangesPos,RangesNeg)
    profileSample <- SummarizedExperiment(profileMat,rowRanges=filteredRanges[match(rownames(profileMat),filteredRanges$giID)])

    ## Set sample name for ChIPprofile object
    
    if(is.null(samplename)){
      if(format %in% c("rlelist","pwm","granges")){
        metadata(profileSample)  <- list(names=c("Sample"))
      }else{
        metadata(profileSample)<- list(names=c(bamFile),AlignedReadsInBam=totalReads)  
      }
    } else{
      metadata(profileSample)<- list(names=samplename,AlignedReadsInBam=totalReads)
    }
    
    ## Pass parameters
    
    paramList <- list("nOfWindows"=nOfWindows,
                      "style"=style,
                      "samplename"=samplename,
                      "nOfWindows"=nOfWindows,
                      "FragmentLength"=FragmentLength,
                      "distanceAround"=distanceAround,
                      "distanceUp"=distanceUp,
                      "distanceDown"=distanceDown,
                      "distanceInRegionStart"=distanceInRegionStart,
                      "distanceInRegionEnd"=distanceInRegionEnd,
                      "distanceOutRegionStart"=distanceOutRegionStart,
                      "distanceOutRegionEnd"=distanceOutRegionEnd,
                      "paired"=paired,
                      "normalize"=normalize,
                      "plotBy"=plotBy,
                      "removeDup"=removeDup,
                      "format"=format,
                      "seqlengths"=seqlengths,
                      "forceFragment"=forceFragment,
                      "method"=method,
                      "genome"=genome,
                      "cutoff"=cutoff,
                      "minFragmentLength"=minFragmentLength,
                      "maxFragmentLength"=maxFragmentLength,
                      "downSample"=downSample
                      )
    return(new("ChIPprofile",profileSample,params=paramList))
  }
  
  ## If style is precentOfRegion. This normalises all regions to the same length.
  
  if(style=="percentOfRegion"){
    
  ##  Spline method average plots
    
    if(method=="spline"){
    
        ### Add option to have different length flanks here?
        
        ## Build extended GRanges of regions and flanks
        
        grWidths <- width(testRangesPos)
        Flanks <- round(grWidths*((distanceAround)/100))      
        RangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-Flanks,end(testRangesPos)+Flanks),strand=Rle("+",length(testRangesPos)),mcols(testRangesPos))     
        grWidths <- width(testRangesNeg)
        Flanks <- round(grWidths*((distanceAround)/100))
        RangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-Flanks,end(testRangesNeg)+Flanks),strand=Rle("+",length(testRangesNeg)),mcols(testRangesNeg))     
        
        ## Initiate empty matrices for counts and GRanges for extracting coverage from rlelist
        matPos <- NULL
        matNeg <- NULL
        testRangesPosNew <- GRanges()
        testRangesNegNew <- GRanges()
        message(paste0("Calculating splines for regions.\nProcessing per contig"))
        
        ## Cycle through contigs and calcualte spline for each view object created from GRanges of interest.
        
        for(c in 1:length(chromosomes)){
          message(paste0("contig: ",i))
          if(any(seqnames(RangesPos)==chromosomes[c])){
            testRangesPosNew <- c(testRangesPosNew,RangesPos[seqnames(RangesPos)==chromosomes[c]
                                                             ])          
          matPos = c(matPos,list(t(viewApply(
                                          Views(genomeCov[names(genomeCov)==chromosomes[c]][[1]],
                                                ranges(RangesPos[seqnames(RangesPos)==chromosomes[c]])
                                          ),
                                            function(x)spline(x,n=(2*(nOfWindows*((distanceAround)/100)))+nOfWindows)$y))))
          }
          
          if(any(seqnames(RangesNeg)==chromosomes[c])){
            testRangesNegNew <- c(testRangesNegNew,RangesNeg[seqnames(RangesNeg)==chromosomes[c]
                                                             ])          
            matNeg = c(matNeg,list(t(viewApply(
              Views(genomeCov[names(genomeCov)==chromosomes[c]][[1]],
                    ranges(RangesNeg[seqnames(RangesNeg)==chromosomes[c]])
                    ),
            function(x)spline(x,n=(2*(nOfWindows*((distanceAround)/100)))+nOfWindows)$y))[,((2*(nOfWindows*((distanceAround)/100)))+nOfWindows):1]))
          }
  
      }
      
      ## Create ChIPprofile object.
      
      message("Creating ChIPprofile")    
      if(!is.null(matPos)){
        matPos <- do.call(rbind,matPos)
      }
      if(!is.null(matNeg)){
        matNeg <- do.call(rbind,matNeg)
      }    
      meansMat <- rbind(matPos,matNeg)
      allRanges <- c(testRangesPosNew,testRangesNegNew)
      rownames(meansMat) <- allRanges$giID
      profileMat <- meansMat[order(rownames(meansMat)),]
      colnames(profileMat) <- c(paste0("Start-",seq(1,(nOfWindows*((distanceAround)/100)))),
                                paste0("Start+",seq(1,nOfWindows)),
                                paste0("End+",seq(1,(nOfWindows*((distanceAround)/100)))))
  
      profileSample <- SummarizedExperiment(profileMat,rowRanges=allRanges[match(rownames(profileMat),allRanges$giID)])
      
      ## Set sample name for ChIPprofile object
      
      if(is.null(samplename)){
        if(format %in% c("rlelist","pwm","granges")){
          metadata(profileSample)  <- list(names=c("Sample"))
        }else{
          metadata(profileSample)<- list(names=c(bamFile),AlignedReadsInBam=totalReads)  
        }
      } else{
        metadata(profileSample)<- list(names=samplename,AlignedReadsInBam=totalReads)
      }
      
      ## Pass parameters
      
      paramList <- list("nOfWindows"=nOfWindows,
                        "style"=style,
                        "samplename"=samplename,
                        "nOfWindows"=nOfWindows,
                        "FragmentLength"=FragmentLength,
                        "distanceAround"=distanceAround,
                        "distanceUp"=distanceUp,
                        "distanceDown"=distanceDown,
                        "distanceInRegionStart"=distanceInRegionStart,
                        "distanceInRegionEnd"=distanceInRegionEnd,
                        "distanceOutRegionStart"=distanceOutRegionStart,
                        "distanceOutRegionEnd"=distanceOutRegionEnd,
                        "paired"=paired,
                        "normalize"=normalize,
                        "plotBy"=plotBy,
                        "removeDup"=removeDup,
                        "format"=format,
                        "seqlengths"=seqlengths,
                        "forceFragment"=forceFragment,
                        "method"=method,
                        "genome"=genome,
                        "cutoff"=cutoff,
                        "minFragmentLength"=minFragmentLength,
                        "maxFragmentLength"=maxFragmentLength,
                        "downSample"=downSample
      )
      return(new("ChIPprofile",profileSample,params=paramList))
    }
      ## Calculate mean coverage within bins across regions.
    if(method=="bin"){
      
      ## 
      meansListNeg <- vector("numeric")
      meansListPos <- vector("numeric")
      
      grListWindowsPos <- GRanges()
      grListWindowsNeg <- GRanges()
      
      ## Create GRanges of windows across regions
      message("Making windows.")
      
      ## Positive regions
      
      if(length(testRangesPos) > 0){
        
        ## Calculate bin lengths
        grWidths <- width(testRangesPos)
        windows <- floor(grWidths%/%nOfWindows)
        extraForWindows <- grWidths%%nOfWindows
        extraForFlankWindows <- grWidths%%(nOfWindows*((distanceAround)/100))
        addToWindow <- 0
        startPos <- start(testRangesPos)-distanceUpStartPos#-(windows/2)
        rem <- rep(0,length(extraForFlankWindows))
        rem2 <- NULL
        
        ## Create bin GRanges for positive 5' flanking regions
        message("Windowing positive 5' flanking ")
        for(i in 1:(nOfWindows*((distanceAround)/100))){
          #message("Window", i,appendLF=F)
          rem2 <- rem+((extraForFlankWindows >= i)+0)
          
          grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(      
            (startPos)+rem+(windows*(i-1)),
            startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))
          rem <- rem2
        }
        
        ## Create bin GRanges for positive regions
        
        startPos <- start(testRangesPos)#-(windows/2)
        rem <- rep(0,length(extraForWindows))
        rem2 <- NULL
        message("Windowing positive regions ")
        for(i in 1:(nOfWindows)){
          #message("Window", i)
          rem2 <- rem+((extraForWindows >= i)+0)
          
          grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(      
            (startPos)+rem+(windows*(i-1)),
            startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))
          rem <- rem2
        }
        
        ## Create bin GRanges for positive 3' flanking regions
        
        rem <- rep(0,length(extraForFlankWindows))
        rem2 <- NULL
        startPos <- end(testRangesPos)#-(windows/2)
        message("Windowing positive 3' flank ")
        for(i in 1:(nOfWindows*((distanceAround)/100))){
          #message("Window", i)
          rem2 <- rem+((extraForFlankWindows >= i)+0)
          
          grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(      
            (startPos)+rem+(windows*(i-1)),
            startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))
          rem <- rem2
        }
        
        ## Order by giID to group windows from gene. Retains secondary order as
        ## created (bin order)
        grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
        
      }
      
       ## Handle negative GRanges as with positive GRanges
      if(length(testRangesNeg) > 0){
        
        grWidths <- width(testRangesNeg)
        windows <- floor(grWidths%/%nOfWindows)
        extraForWindows <- grWidths%%nOfWindows
        extraForFlankWindows <- grWidths%%(nOfWindows*((distanceAround)/100))
        addToWindow <- 0
        startPos <- start(testRangesNeg)-distanceDownEndNeg
        rem <- rep(0,length(extraForFlankWindows))
        rem2 <- NULL
        message("Windowing negative 5' flank ")
        for(i in 1:(nOfWindows*((distanceAround)/100))){
          #message("Window", i)
          rem2 <- rem+((extraForFlankWindows >= i)+0)
          
          grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(      
            (startPos)+rem+(windows*(i-1)),
            startPos+(windows*i)-1+rem2),giID=testRangesNeg$giID))
          rem <- rem2
        }
        
        startPos <- start(testRangesNeg)#-(windows/2)
        rem <- rep(0,length(extraForWindows))
        rem2 <- NULL
        message("Windowing negative regions ")
        
        for(i in 1:(nOfWindows)){
          #message("Window", i)
          rem2 <- rem+((extraForWindows >= i)+0)
          
          grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(      
            (startPos)+rem+(windows*(i-1)),
            startPos+(windows*i)-1+rem2),giID=testRangesNeg$giID))
          rem <- rem2
        }
        rem <- rep(0,length(extraForFlankWindows))
        rem2 <- NULL
        startPos <- end(testRangesNeg)#-(windows/2)
        message("Windowing negative 3' flank ")
        
        for(i in 1:(nOfWindows*((distanceAround)/100))){
          #message("Window", i)
          rem2 <- rem+((extraForFlankWindows >= i)+0)
          
          grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(      
            (startPos)+rem+(windows*(i-1)),
            startPos+(windows*i)-1+rem2),giID=testRangesNeg$giID))
          rem <- rem2
        }
        
        grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
      
      }
      grListWindows <- list(grListWindowsPos,grListWindowsNeg)
      message("..done\n")
      
      ## Cycle through contigs to extract scores from rlelist per contig

      message(paste0("Calculating bin scores for regions.\nProcessing per contig"))
      
      
      for(c in 1:length(chromosomes)){
 
        message(paste0("contig: ",c))
        
        message("Processing inner region windows in ",chromosomes[c])
        covPerPeakPos <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[[1]][seqnames(grListWindows[[1]]) == chromosomes[c]]))
        doubleTempPos <- viewMeans(covPerPeakPos)
        names(doubleTempPos) <- as.vector(grListWindows[[1]][seqnames(grListWindows[[1]]) == chromosomes[c]]$giID)
        meansListPos <- c(meansListPos,doubleTempPos)
        covPerPeakNeg <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[[2]][seqnames(grListWindows[[2]]) == chromosomes[c]]))
        doubleTempNeg <- viewMeans(covPerPeakNeg)
        names(doubleTempNeg) <- as.vector(grListWindows[[2]][seqnames(grListWindows[[2]]) == chromosomes[c]]$giID)
        meansListNeg <- c(meansListNeg,doubleTempNeg)
        message("..done")
        message("Processing flanking windows in ",chromosomes[c])      
        
        tempstartRegionRangesPosMat <- NULL
        tempendRegionRangesPosMat <- NULL  
        tempstartRegionRangesNegMat <- NULL
        tempendRegionRangesNegMat <- NULL    
      }
      
      ## Create matrices for mean bin coverage
      
      meansPos <- matrix(meansListPos,
                             ncol=((nOfWindows*((distanceAround)/100))*2)+nOfWindows
                             ,byrow=TRUE)
      if(nrow(meansPos) > 0){
      rownames(meansPos) <- matrix(names(meansListPos),ncol=((nOfWindows*((distanceAround)/100))*2)+nOfWindows
                                       ,byrow=TRUE)[,1]
      }
      meansNeg <- matrix(meansListNeg,
                         ncol=((nOfWindows*((distanceAround)/100))*2)+nOfWindows
                         ,byrow=TRUE)[,(((nOfWindows*((distanceAround)/100))*2)+nOfWindows):1]
      if(nrow(meansNeg) > 0){
        
      rownames(meansNeg) <- matrix(names(meansListNeg),ncol=((nOfWindows*((distanceAround)/100))*2)+nOfWindows
                                   ,byrow=TRUE)[,1]
      }
      meansMat <- rbind(meansPos,meansNeg)
      profileMat <- meansMat[order(rownames(meansMat)),]
    }
    
    ### Create ChIPprofile object.
    
    colnames(profileMat) <- c(paste0("Start-",seq(1,(nOfWindows*((distanceAround)/100)))),
                              paste0("Start+",seq(1,nOfWindows)),
                              paste0("End+",seq(1,(nOfWindows*((distanceAround)/100)))))
    filteredRanges <- c(testRangesPos,testRangesNeg)

    profileSample <- SummarizedExperiment(profileMat,rowRanges=filteredRanges[match(rownames(profileMat),filteredRanges$giID)])

    
    if(is.null(samplename)){
      if(format %in% c("rlelist","pwm","granges")){
        metadata(profileSample)  <- list(names=c("Sample"))
      }else{
        metadata(profileSample)<- list(names=c(bamFile),AlignedReadsInBam=totalReads)  
      }
    }else{
      metadata(profileSample)<- list(names=samplename,AlignedReadsInBam=totalReads)
    }
    
    ## Pass parameters
    
    paramList <- list("nOfWindows"=nOfWindows,
                      "style"=style,
                      "samplename"=samplename,
                      "nOfWindows"=nOfWindows,
                      "FragmentLength"=FragmentLength,
                      "distanceAround"=distanceAround,
                      "distanceUp"=distanceUp,
                      "distanceDown"=distanceDown,
                      "distanceInRegionStart"=distanceInRegionStart,
                      "distanceInRegionEnd"=distanceInRegionEnd,
                      "distanceOutRegionStart"=distanceOutRegionStart,
                      "distanceOutRegionEnd"=distanceOutRegionEnd,
                      "paired"=paired,
                      "normalize"=normalize,
                      "plotBy"=plotBy,
                      "removeDup"=removeDup,
                      "format"=format,
                      "seqlengths"=seqlengths,
                      "forceFragment"=forceFragment,
                      "method"=method,
                      "genome"=genome,
                      "cutoff"=cutoff,
                      "minFragmentLength"=minFragmentLength,
                      "maxFragmentLength"=maxFragmentLength,
                      "downSample"=downSample
    )
    return(new("ChIPprofile",profileSample,params=paramList))
  } 

  ## Run when style is Region. This style creates a hybrid plot where the edges are presented
  ## as per base pair but the remaining region is binned to common normalised length.
  
  if(style=="region"){
    
    message("Defining flanks of regions..",appendLF=FALSE)
    
    ## Create GRanges for flanking regions
    startRegionRangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)-distanceOutRegionStart,start(testRangesPos)+distanceInRegionStart),strand=Rle("+",length(testRangesPos)),mcols(testRangesPos))
    endRegionRangesPos <- GRanges(seqnames(testRangesPos),IRanges(end(testRangesPos)-distanceInRegionEnd,end(testRangesPos)+distanceOutRegionEnd),strand=Rle("+",length(testRangesPos)),mcols(testRangesPos))
    startRegionRangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(end(testRangesNeg)-distanceInRegionStart,end(testRangesNeg)+distanceOutRegionStart),strand=Rle("+",length(testRangesNeg)),mcols(testRangesNeg))
    endRegionRangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)-distanceOutRegionEnd,start(testRangesNeg)+distanceInRegionEnd),strand=Rle("+",length(testRangesNeg)),mcols(testRangesNeg))
    
    testRangesPos <- GRanges(seqnames(testRangesPos),IRanges(start(testRangesPos)+distanceInRegionStart,end(testRangesPos)-distanceInRegionEnd),strand=Rle("+",length(testRangesPos)),mcols(testRangesPos))
    testRangesNeg <- GRanges(seqnames(testRangesNeg),IRanges(start(testRangesNeg)+distanceInRegionEnd,end(testRangesNeg)-distanceInRegionStart),strand=Rle("+",length(testRangesNeg)),mcols(testRangesNeg))     
    
    message("...Done")
    
    ## Create windows for regions (minus flanking)
    ## Same method as for percentOfRegion
    
    meansList <- vector("numeric")
    grListWindowsPos <- GRanges()
    grListWindowsNeg <- GRanges()
    message("Making windows")
    
    
    
    if(length(testRangesPos) > 0){
      grWidths <- width(testRangesPos)  
      windows <- floor(grWidths%/%nOfWindows)
      
      ########################
#      extraLastWindow <- grWidths%%nOfWindows
#      addToWindow <- 0
      
      ### This is done differently to percentOfRegion style.?
      
      ## Create GRanges windows of binned regions
      
#      for(i in 1:nOfWindows){
        
#        if(i == nOfWindows){
#          addToWindow <- extraLastWindow 
#        }
#        grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(      
#          start(testRangesPos)+(windows*(i-1)),
#          start(testRangesPos)+(windows*i)-1+addToWindow),giID=testRangesPos$giID))
#      }
#    }
#    grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
    

    ###########
  
    
    extraForWindows <- grWidths%%nOfWindows
    addToWindow <- 0
    startPos <- start(testRangesPos)#-(windows/2)
    rem <- rep(0,length(extraForWindows))
    rem2 <- NULL
    message("Windowing positive regions ")
    for(i in 1:(nOfWindows)){
      #message("Window", i)
      rem2 <- rem+((extraForWindows >= i)+0)
      
      grListWindowsPos <- c(grListWindowsPos,GRanges(seqnames(testRangesPos),IRanges(      
        (startPos)+rem+(windows*(i-1)),
        startPos+(windows*i)-1+rem2),giID=testRangesPos$giID))
      rem <- rem2
    }
    grListWindowsPos <- grListWindowsPos[order(grListWindowsPos$giID)]
    
    
    ############
    
    
    if(length(testRangesNeg) > 0){
      grWidths <- width(testRangesNeg)
      windows <- floor(grWidths%/%nOfWindows)
#      extraLastWindow <- grWidths%%nOfWindows
#      addToWindow <- 0
# 
#       for(i in 1:nOfWindows){
#         
#         if(i == nOfWindows){
#           addToWindow <- extraLastWindow 
#         }
#         grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(     
#           end(testRangesNeg)-(windows*i)+1-addToWindow,
#           end(testRangesNeg)-(windows*(i-1))),giID=testRangesNeg$giID))
#       }
#       grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
###########

###########


    extraForWindows <- grWidths%%nOfWindows
    addToWindow <- 0
    startNeg <- start(testRangesNeg)#-(windows/2)
    rem <- rep(0,length(extraForWindows))
    rem2 <- NULL
    message("Windowing negative regions ")
    for(i in 1:(nOfWindows)){
      #message("Window", i)
      rem2 <- rem+((extraForWindows >= i)+0)
      
      grListWindowsNeg <- c(grListWindowsNeg,GRanges(seqnames(testRangesNeg),IRanges(      
        (startNeg)+rem+(windows*(i-1)),
        startNeg+(windows*i)-1+rem2),giID=testRangesNeg$giID))
      rem <- rem2
    }
    grListWindowsNeg <- grListWindowsNeg[order(grListWindowsNeg$giID)]
    

############

############


    
    }
    grListWindows <- c(grListWindowsPos,grListWindowsNeg)

    message(paste0("Calculating bin scores for regions and per base pair for flanks.                   
                   \nProcessing per contig"))
  
    for(c in 1:length(chromosomes)){
      
      message(paste0("contig: ",c))
      
      ## Extract mean coverage of Views from Granges of regions
      
     # message("Processing inner region windows in ",chromosomes[c])
      covPerPeak <- Views(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]],ranges(grListWindows[seqnames(grListWindows) == chromosomes[c]]))
      doubleTemp <- viewMeans(covPerPeak)
      names(doubleTemp) <- as.vector(grListWindows[seqnames(grListWindows) == chromosomes[c]]$giID)
      meansList <- c(meansList,doubleTemp)
      #message("..done")
      #message("Processing flanking windows in ",chromosomes[c])      
      
      ## Perform a point style retrieval of signal from around edges.
     
      tempstartRegionRangesPosMat <- NULL
      tempendRegionRangesPosMat <- NULL  
      tempstartRegionRangesNegMat <- NULL
      tempendRegionRangesNegMat <- NULL    
      
      if(length(startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]]) > 0){
        tempstartRegionRangesPosMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]])]),ncol=mean(width(startRegionRangesPos)),byrow=TRUE)
        rownames(tempstartRegionRangesPosMat) <- startRegionRangesPos[seqnames(startRegionRangesPos) %in% chromosomes[c]]$giID
      }
      
      if(length(endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]]) > 0){
        tempendRegionRangesPosMat <- matrix(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]])]),ncol=mean(width(endRegionRangesPos)),byrow=TRUE)
        rownames(tempendRegionRangesPosMat) <- endRegionRangesPos[seqnames(endRegionRangesPos) %in% chromosomes[c]]$giID
      }
      if(length(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]]) > 0){
        tempstartRegionRangesNegMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]])])),ncol=mean(width(startRegionRangesNeg)),byrow=TRUE)
        rownames(tempstartRegionRangesNegMat) <- rev(startRegionRangesNeg[seqnames(startRegionRangesNeg) %in% chromosomes[c]]$giID)
      }
      if(length(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]]) > 0){
        tempendRegionRangesNegMat <- matrix(rev(as.vector(genomeCov[[which(names(genomeCov) %in% chromosomes[c])]][ranges(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]])])),ncol=mean(width(endRegionRangesNeg)),byrow=TRUE)
        rownames(tempendRegionRangesNegMat) <- rev(endRegionRangesNeg[seqnames(endRegionRangesNeg) %in% chromosomes[c]]$giID)
      }
      
      posRegionStartMat <- rbind(posRegionStartMat,tempstartRegionRangesPosMat)
      posRegionEndMat <- rbind(posRegionEndMat,tempendRegionRangesPosMat)
      negRegionStartMat <- rbind(negRegionStartMat,tempstartRegionRangesNegMat)
      negRegionEndMat <- rbind(negRegionEndMat,tempendRegionRangesNegMat)
      tempstartRegionRangesPosMat <- NULL
      tempendRegionRangesPosMat <- NULL  
      tempstartRegionRangesNegMat <- NULL
      tempendRegionRangesNegMat <- NULL           
      message("..done")
    }

    ## Create ChIP profile
    message("Creating ChIPprofile")

    AllRegionStart <- rbind(posRegionStartMat,negRegionStartMat)
    AllRegionEnd <- rbind(posRegionEndMat,negRegionEndMat)
    meansMat <- matrix(meansList,ncol=nOfWindows,byrow=TRUE)
    rownames(meansMat) <- matrix(names(meansList),ncol=nOfWindows,byrow=TRUE)[,1]
    start <- cbind(seq(1,length(colMeans(AllRegionStart))),colMeans(AllRegionStart))
    mid <- cbind(max(start[,1])+seq(1,length(colMeans(meansMat)))*nOfWindows,colMeans(meansMat))
    end <- cbind(max(mid[,1])+seq(1,length(colMeans(AllRegionEnd))),colMeans(AllRegionEnd))
    profileMat <- cbind(AllRegionStart[order(rownames(AllRegionStart)),],
                        meansMat[order(rownames(meansMat)),],
                        AllRegionEnd[order(rownames(AllRegionEnd)),])
    colnames(profileMat) <- c(paste0("Region_Start",seq(0-distanceOutRegionStart,-1)),"Region_Start",paste0("Region_Start",seq(1,distanceInRegionStart)),
    paste0(seq(1,nOfWindows),"%_ofRegion"),
    paste0("Region_End",seq(0-distanceInRegionEnd,-1)),"Region_End",paste0("Region_End",seq(1,distanceOutRegionEnd))
    )
    filteredRanges <- c(testRangesPos,testRangesNeg)
    profileSample <- SummarizedExperiment(profileMat,rowRanges=filteredRanges[match(rownames(profileMat),filteredRanges$giID)])
    print(format)

    if(is.null(samplename)){
      if(format %in% c("rlelist","pwm","granges")){
        metadata(profileSample)  <- list(names=c("Sample"))
      }else{
        metadata(profileSample)<- list(names=c(bamFile),AlignedReadsInBam=totalReads)  
      }
    }else{
      metadata(profileSample)<- list(names=samplename,AlignedReadsInBam=totalReads)
    }
    
    ## Pass parameters
    
    paramList <- list("nOfWindows"=nOfWindows,
                      "style"=style,
                      "samplename"=samplename,
                      "nOfWindows"=nOfWindows,
                      "FragmentLength"=FragmentLength,
                      "distanceAround"=distanceAround,
                      "distanceUp"=distanceUp,
                      "distanceDown"=distanceDown,
                      "distanceInRegionStart"=distanceInRegionStart,
                      "distanceInRegionEnd"=distanceInRegionEnd,
                      "distanceOutRegionStart"=distanceOutRegionStart,
                      "distanceOutRegionEnd"=distanceOutRegionEnd,
                      "paired"=paired,
                      "normalize"=normalize,
                      "plotBy"=plotBy,
                      "removeDup"=removeDup,
                      "format"=format,
                      "seqlengths"=seqlengths,
                      "forceFragment"=forceFragment,
                      "method"=method,
                      "genome"=genome,
                      "cutoff"=cutoff,
                      "minFragmentLength"=minFragmentLength,
                      "maxFragmentLength"=maxFragmentLength,
                      "downSample"=downSample
    )
    return(new("ChIPprofile",profileSample,params=paramList))
        
    } 

}
}

