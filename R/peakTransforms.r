#' Plot coverage of points or regions.
#' 
#' @rdname findconsensusRegions
#' @param testRanges Named character vector of region locations 
#' @param bamFiles Named character vector of bamFile locations
#' @param method Method to select reproducible summits to merge.
#' @param summit Only mean avaialble
#' @param resizepeak Only asw available
#' @param overlap Type of overlap to consider for finding consensus sites
#' @param fragmentLength Predicted fragment length. Set to NULL to auto-calculate
#' @param NonPrimaryPeaks A list of parameters to deal with non primary peaks in consensus regions.
#' @return Consensus A GRanges object of consensus regions with consensus summits. 
#' @export
findconsensusRegions <- function(testRanges,bamFiles=NULL,method="majority",summit="mean",resizepeak="asw",overlap="any",fragmentLength=NULL,
                                 NonPrimaryPeaks=list(withinsample="drop",betweensample="mean")){

  testRanges <- GRangesList(
    bplapply(
      testRanges,
      GetGRanges)
    )
  
  ans <- lapply(
    names(bamFiles),
    function(x) summitPipeline(
      unlist(bamFiles[x]),
      unlist(testRanges[x]),
      fragmentLength=NULL,readlength=36
      )
  )
  consensusRanges <- runConsensusRegions(
    testRanges,
    method="majority",
    overlap="any"
    )
  if(unlist(NonPrimaryPeaks["withinsample"])=="drop"){
    consensensusAns <- bplapply(ans,function(x)
      dropNonPrimary(x,consensusRanges)
    )
    ansSummits <- do.call(cbind,bplapply(ans,function(x)
      extractSummits(x,consensusRanges)
    ))
    ansSummitScores <- do.call(cbind,bplapply(ans,function(x)
      extractScores(x,consensusRanges)
    )) 
    if(unlist(NonPrimaryPeaks["betweensample"])=="mean"){
      meanSummits <- rowMeans(ansSummits,na.rm=TRUE)
    }
    if(unlist(NonPrimaryPeaks["betweensample"])=="weightedmean"){
      #meanSummits <- rowMeans(ansSummits,na.omit=TRUE)
      meanSummits <- round(sapply(seq(1,nrow(ansSummits)),function(x)weighted.mean(ansSummits[x,],ansSummitScores[x,])))
    }
  }  
  start(consensusRanges) <- end(consensusRanges) <- meanSummits
  return(consensusRanges)  
}

extractSummits <- function(x,consensusRanges){
  summits <- vector("numeric",length=length(consensusRanges))
  overMat <- as.matrix(findOverlaps(consensusRanges,x))
  summits[overMat[,1]] <- as.vector(start(x[overMat[,2]]))
  return(summits)
}
extractScores <- function(x,consensusRanges,score="summitScores"){
  scores <- vector("numeric",length=length(consensusRanges))
  overMat <- as.matrix(findOverlaps(consensusRanges,x))
  scores[overMat[,1]] <- as.vector(mcols(x[overMat[,2]])[,"summitScores"])
  return(scores)
}

dropNonPrimary <- function(x,consensusRanges,id="mcols.ID",score="summitScores"){
  mat <- as.matrix(findOverlaps(consensusRanges,x))
  tempConsensusRanges <- consensusRanges[mat[,1],]
  mcols(tempConsensusRanges) <- mcols(x[mat[,2]])[,c(id,score)]
  tempConsensusRanges <- tempConsensusRanges[order(mcols(tempConsensusRanges)[,score],decreasing=TRUE),]  
  primaryIDs <- mcols(tempConsensusRanges[match(unique(tempConsensusRanges[,id]),tempConsensusRanges[,id])])[,id]
  x <- x[mcols(x)[,id] %in% primaryIDs]
  return(x)
}
#' Returns summits and summmit scores after optional fragment length prediction and read extension
#' @rdname findconsensusRegions
#' @param peakfile GRanges of genomic intervals to summit. 
#' @param reads Character vector of bamFile location or GAlignments object
#' @param readlength Read length of alignments.
#' @return Summits A GRanges object of summits and summit scores.
#' @export
summitPipeline <- function(reads,peakfile,fragmentLength,readlength){
  message("Reading in peaks..",appendLF=FALSE)
  testRanges <- GetGRanges(peakfile)
  message("done")  
  if(class(reads) == "GAlignments"){
    message("Alignments loaded")
    ChrLengths <- seqlengths(reads)
  }
  if(class(reads) == "character"){
    message("Reading in alignments..",appendLF=FALSE)
    ChrLengths <- scanBamHeader(reads)[[1]]$targets
    reads <- readGAlignmentsFromBam(reads)
    message("done") 
  }
  message("Calculating fragmentlength..",appendLF=FALSE)
  ccscores <- getShifts(reads,ChrLengths,shiftWindowStart=1,shiftWindowEnd=400)
  fragmentLength <- getFragmentLength(ccscores,readlength)
  message("done")
  message("Extending reads..",appendLF=FALSE)  
  reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
  message("done")
  message("Finding summit locations..",appendLF=FALSE)
  peaks <- runFindSummit(testRanges,reads,fragmentLength=NULL)
  message("done")
  message("Scoring summits..",appendLF=FALSE)
  peaks <- getSummitScore(reads,peaks,fragmentLength=NULL,score="height")
  message("done")
  return(peaks)
}

runConsensusRegions <- function(testRanges,method="majority",overlap="any"){
    if(class(testRanges) == "GRangesList" & length(testRanges) > 1){
      
      reduced <- reduce(unlist(testRanges))
      consensusIDs <- paste0("consensus_",seq(1,length(reduced)))
      mcols(reduced) <- 
      do.call(cbind,lapply(testRanges,function(x)(reduced %over% x)+0))
      if(method=="majority"){
        reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) > length(testRanges)/2,]
      }
      if(method=="none"){
        reducedConsensus <- reduced
      }
      if(is.numeric(method)){
        reducedConsensus <- reduced[rowSums(as.data.frame(mcols(reduced))) > method,]
      }
    consensusIDs <- paste0("consensus_",seq(1,length(reducedConsensus)))
    mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)),consensusIDs)
    return(reducedConsensus)
    
  }
}
runGetShifts <- function(reads,ChrLengths,ChrOfInterestshift,shiftWindowStart=1,shiftWindowEnd=400){

reads <- reads
ChrLengths <- seqlengths(reads)
PosCoverage <- coverage(IRanges(start(reads[strand(reads)=="+"]),start(reads[strand(reads)=="+"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
NegCoverage <- coverage(IRanges(end(reads[strand(reads)=="-"]),end(reads[strand(reads)=="-"])),width=ChrLengths[names(ChrLengths) %in% ChrOfInterestshift])
message("Calculating shift for ",ChrOfInterestshift,"\n")
ShiftsTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),PosCoverage,NegCoverage,RleSumAny, verbose = TRUE)         
return(ShiftsTemp)
}
getShifts <- function(reads,ChrLengths,
                      shiftWindowStart=1,shiftWindowEnd=400){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }  
shiftMat <- do.call(cbind,bplapply(names(ChrLengths),function(x)
runGetShifts(reads[seqnames(reads) %in% x],ChrLengths,x,
          shiftWindowStart=1,shiftWindowEnd=400)))
cc_scores <- (rowSums(shiftMat)[1]-rowSums(shiftMat))/rowSums(shiftMat)[1]
return(cc_scores)
}

getFragmentLength <- function(x,readLength){
  #peaks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0)+1
  MaxShift <- which.max(runMean(x[-seq(1,(2*readLength))],10))+2*readLength
  
}
  
runFindSummit <- function(testRanges,reads,fragmentLength=NULL){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }
  if(!is.null(fragmentLength)){
    message("Extending reads to fragmentlength of ",fragmentLength,"..",appendLF=FALSE)
    reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
    message("done")
  }
  test <- do.call(c,
                  bplapply(
    unique(seqnames(reads))[unique(seqnames(reads)) %in% unique(seqnames(testRanges))],
    function(x) 
    findCovMaxPos(reads[seqnames(testRanges) %in% x],testRanges[seqnames(testRanges) %in% x],seqlengths(reads)[names(seqlengths(reads)) %in% x],fragmentLength)
    )
  )
  return(test)                                        
}

getSummitScore <- function(reads,summits,fragmentLength=NULL,score="height"){
  if(is.character(reads)){
    reads <- readGAlignmentsFromBam(reads)
  }
  if(!is.null(fragmentLength)){
    message("Extending reads to fragmentlength of ",fragmentLength,"..",appendLF=FALSE)
    reads <- resize(as(reads,"GRanges"),fragmentLength,"start")
    message("done")
  }
  test <- do.call(c,
                  bplapply(
                    as.vector(unique(seqnames(reads))[unique(seqnames(reads)) %in% unique(seqnames(summits))]),
                    function(x) 
                      runGetSummitScore(reads[seqnames(reads) %in% x],summits[seqnames(summits) %in% x],seqlengths(reads)[names(seqlengths(reads)) %in% x])
                  )
  )
  return(test)                                         
}
runGetSummitScore <- function(reads,summits,ChrOfInterestshift,FragmentLength=150,score="height"){
    

  cov <- coverage(reads)
    if(score=="height"){
      summitScores <- as.vector(unlist(cov[summits],use.names=FALSE))
    }
  mcols(summits) <- cbind(as.data.frame(mcols(summits)),summitScores)
    return(summits)
}

RleSumAny <- function (e1, e2)
{
  len <- length(e1)
  stopifnot(len == length(e2))
  x1 <- runValue(e1); s1 <- cumsum(runLength(e1))
  x2 <- runValue(e2); s2 <- cumsum(runLength(e2))
  .Call("rle_sum_any",
        as.integer(x1), as.integer(s1),
        as.integer(x2), as.integer(s2),
        as.integer(len),
        PACKAGE = "chipseq")
}
library(BiocParallel)
library(GenomicAlignments)


#' Create GRangeslist from all combinations of GRanges
#'
#' @param testRanges A named list of GRanges or a named GRangesList 
#' @return groupedGRanges A named GRangesList object.
#' @examples 
#' data(ik_Example)
#'  gts <- groupByOverlaps(ik_Example)
#' @export
groupByOverlaps <- function(testRanges){
  
  testRanges <- GRangesList(
    lapply(
      testRanges,
      GetGRanges)
  )
  allRegionsReduced <- reduce(
    unlist(testRanges)
  )
  
  overlapMat <- do.call(cbind,
                        lapply(1:length(testRanges),
                               function(x) ifelse((allRegionsReduced %over% testRanges[[x]]),names(testRanges)[x],""))
                        )
  
  colnames(overlapMat) <- names(testRanges)
  mcols(allRegionsReduced)$grangesGroups <- as.factor(gsub("--|^-|-$","",
                                                 apply(overlapMat, 1 , paste , collapse = "-" )
                                                ))
  
  groupedGRangesList <- lapply(levels(allRegionsReduced$grangesGroups),
         function(x)allRegionsReduced[allRegionsReduced$grangesGroups %in% x])
  names(groupedGRangesList) <- levels(allRegionsReduced$grangesGroups)
  return(groupedGRangesList)
}

#' Set strand by overlapping or nearest anchor GRanges
#' @rdname orientBy
#' @param testRanges The GRanges object to anchor. 
#' @param anchorRanges A GRanges object by which to anchor strand orientation. 
#' @return newRanges A GRanges object.
#' @examples
#' data(ik_Example)
#' strand(ik_Example[[1]]) <- "+"
#' anchoredGRanges <- orientBy(ik_Example[[2]],ik_Example[[1]]) 
#' @export
orientBy <- function(testRanges,anchorRanges){
  distIndex <- distanceToNearest(testRanges,anchorRanges)
  anchorRangesFilt <- anchorRanges[subjectHits(distIndex)]
#   widths <- width(pintersect(testRanges,anchorRangesFilt,resolve.empty="max.start"))
#   testRanges$distance <- mcols(distIndex)$distance
#   testRanges$overlapsize <- widths
#   anchorRangesFilt <- anchorRangesFilt[order(mcols(distIndex)$distance,widths),]
  strand(testRanges) <- strand(anchorRangesFilt)
  return(testRanges)
}


GetGRanges <- function(LoadFile,AllChr=NULL,ChrOfInterest=NULL,simple=FALSE,sepr="\t",simplify=FALSE){
  #    require(Rsamtools)
  #    require(GenomicRanges)
  
  if(class(LoadFile) == "GRanges"){
    RegionRanges <- LoadFile
    if(simplify){
      RegionRanges <- GRanges(seqnames(RegionRanges),ranges(RegionRanges))
    }
  }else{
    if(class(LoadFile) == "character"){
      RangesTable <- read.delim(LoadFile,sep=sepr,header=TRUE,comment.char="#")
    }else if(class(LoadFile) == "matrix"){
      RangesTable <- as.data.frame(LoadFile)
    } else{
      RangesTable <- as.data.frame(LoadFile)
    }
    Chromosomes <- as.vector(RangesTable[,1])
    Start <- as.numeric(as.vector(RangesTable[,2]))
    End <- as.numeric(as.vector(RangesTable[,3]))
    RegionRanges <- GRanges(seqnames=Chromosomes,ranges=IRanges(start=Start,end=End))
    if(simple == FALSE){
      if(ncol(RangesTable) > 4){
        ID <- as.vector(RangesTable[,4])
        Score <- as.vector(RangesTable[,5])
        if(ncol(RangesTable) > 6){
          Strand <- rep("*",nrow(RangesTable))
          RemainderColumn <- as.data.frame(RangesTable[,-c(1:6)])
          mcols(RegionRanges) <- cbind(ID,Score,Strand,RemainderColumn)
        }else{
          mcols(RegionRanges) <- cbind(ID,Score)
        }
      }
    }
  }
  if(!is.null(AllChr)){ 
    RegionRanges <- RegionRanges[seqnames(RegionRanges) %in% AllChr]    
    seqlevels(RegionRanges,force=TRUE) <- AllChr
  }
  if(!is.null(ChrOfInterest)){      
    RegionRanges <- RegionRanges[seqnames(RegionRanges) == ChrOfInterest]      
  }
  
  return(RegionRanges)
}

findCovMaxPos <- function(reads,bedRanges,ChrOfInterest,FragmentLength){
  #    require(GenomicRanges)
  #    require(Rsamtools)
  
  cat("done\n")
  cat("Calculating coverage\n")
  MaxRanges <- GRanges()
  if(length(reads) > 0){
    seqlengths(reads)[names(ChrOfInterest)] <- ChrOfInterest
    AllCov <- coverage(reads) 
    cat("Calculating Summits on ",names(ChrOfInterest)," ..")
    covPerPeak <- Views(AllCov[[which(names(AllCov) %in% names(ChrOfInterest))]],ranges(bedRanges[seqnames(bedRanges) == names(ChrOfInterest)]))
    meanSummitLocations <- viewApply(covPerPeak,function(x)round(mean(which(x==max(x)))))
    Maxes <- (start(bedRanges)+meanSummitLocations)-1
    if(any(is.na(Maxes))){ 
      NoSummitRanges <- bedRanges[is.na(Maxes)]
      Maxes[is.na(Maxes)]  <- (start((ranges(NoSummitRanges[seqnames(NoSummitRanges) == names(ChrOfInterest)])))+end((ranges(NoSummitRanges[seqnames(NoSummitRanges) == names(ChrOfInterest)]))))/2
    }
    MaxRanges <- GRanges(seqnames(bedRanges[seqnames(bedRanges) == names(ChrOfInterest)]),IRanges(start=Maxes,end=Maxes),mcols=mcols(bedRanges[seqnames(bedRanges) == names(ChrOfInterest)]))
    #revAllCov <- rev(coverage(reads))
    #revAllCov <- runmean(revAllCov[names(revAllCov) %in% ChrOfInterest],20)
    #cat("Calculating reverse Summits on ",ChrOfInterest," ..")
    #revMaxes <- which.max(Views(revAllCov[[which(names(revAllCov) %in% ChrOfInterest)]],ranges(bedRanges[seqnames(bedRanges) == ChrOfInterest])))
    #if(any(is.na(revMaxes))){ 
    #  revNoSummitRanges <- bedRanges[is.na(revMaxes)]
    #  revMaxes[is.na(revMaxes)]  <- (start((ranges(revNoSummitRanges[seqnames(revNoSummitRanges) == ChrOfInterest])))+end((ranges(revNoSummitRanges[seqnames(revNoSummitRanges) == ChrOfInterest]))))/2
    #}
    #revMaxRanges <- GRanges(seqnames(bedRanges[seqnames(bedRanges) == ChrOfInterest]),IRanges(start=Maxes,end=Maxes),mcols=mcols(bedRanges[seqnames(bedRanges) == ChrOfInterest]))
    #meanMaxes <- rowMeans(cbind(Maxes,revMaxes))
    #meanMaxRanges <- GRanges(seqnames(bedRanges[seqnames(bedRanges) == ChrOfInterest]),IRanges(start=meanMaxes,end=meanMaxes),mcols=mcols(bedRanges[seqnames(bedRanges) == ChrOfInterest]))
    #cat(".done\n")
  }
  #return(meanMaxRanges)
  return(MaxRanges)
}

runMean = function(x, k, alg=c("C", "R", "fast", "exact"),
                   endrule=c("mean", "NA", "trim", "keep", "constant", "func"),
                   align = c("center", "left", "right"))
{
  alg     = match.arg(alg)
  endrule = match.arg(endrule)
  align   = match.arg(align)
  dimx = dim(x) # Capture dimension of input array - to be used for formating y
  x = as.vector(x)
  n = length(x)
  if (k<=1) return (x)
  if (k >n) k = n
  k2 = k%/%2
    y = double(n)
    k1 = k-k2-1
    y = c( sum(x[1:k]), diff(x,k) ); # find the first sum and the differences from it
    y = cumsum(y)/k                  # apply precomputed differences
    y = c(rep(0,k1), y, rep(0,k2))   # make y the same length as x
    if (endrule=="mean") endrule="func"
  y = EndRule(x, y, k, dimx, endrule, align, mean, na.rm=TRUE)
  return(y)
}

EndRule = function(x, y, k, dimx,
                   endrule=c("NA", "trim", "keep", "constant", "func"),
                   align = c("center", "left", "right"), Func, ...)
{
  # Function which postprocess results of running windows functions and cast
  # them in to specified format. On input y is equivalent to
  #   y = runFUNC(as.vector(x), k, endrule="func", align="center")
  
  # === Step 1: inspects inputs and unify format ===
  align   = match.arg(align)
  k = as.integer(k)
  k2 = k%/%2
  if (k2<1) k2 = 1
  yIsVec = is.null(dimx) # original x was a vector -> returned y will be a vector
  if (yIsVec) dimx=c(length(y),1) # x & y will become 2D arrays
  dim(x) <- dimx
  dim(y) <- dimx
  n = nrow(x)
  m = ncol(x)
  if (k>n) k2 = (n-1)%/%2
  k1 = k-k2-1
  if (align=="center" && k==2) align='right'
  
  # === Step 2: Apply different endrules ===
  if (endrule=="trim") {
    y = y[(k1+1):(n-k2),] # change y dimensions
  } else if (align=="center") {
    idx1 = 1:k1
    idx2 = (n-k2+1):n
    # endrule calculation in R will be skipped for most common case when endrule
    # is default and array was a vector not a matrix
    if (endrule=="NA") {
      y[idx1,] = NA
      y[idx2,] = NA
    } else if (endrule=="keep") {
      y[idx1,] = x[idx1,]
      y[idx2,] = x[idx2,]
    } else if (endrule=="constant") {
      y[idx1,] = y[k1+1+integer(m),]
      y[idx2,] = y[n-k2+integer(m),]
    } else if (endrule=="func" || !yIsVec) {
      for (j in 1:m) {
        for (i in idx1) y[i,j] = Func(x[1:(i+k2),j], ...)
        for (i in idx2) y[i,j] = Func(x[(i-k1):n,j], ...)
      }
    }
  } else if (align=="left") {
    y[1:(n-k1),] = y[(k1+1):n,]
    idx = (n-k+2):n
    if (endrule=="NA") {
      y[idx,] = NA
    } else if (endrule=="keep") {
      y[idx,] = x[idx,]
    } else if (endrule=="constant") {
      y[idx,] = y[n-k+integer(m)+1,]
    } else {
      for (j in 1:m) for (i in idx) y[i,j] = Func(x[i:n,j], ...)
    }
  } else if (align=="right") {
    y[(k2+1):n,] = y[1:(n-k2),]
    idx = 1:(k-1)
    if (endrule=="NA") {
      y[idx,] = NA
    } else if (endrule=="keep") {
      y[idx,] = x[idx,]
    } else if (endrule=="constant") {
      y[idx,] = y[k+integer(m),]
    } else {
      for (j in 1:m) for (i in idx) y[i,j] = Func(x[1:i,j], ...)
    }
  }
  
  # === Step 4: final casting and return results ===
  if (yIsVec) y = as.vector(y);
  return(y)
}

