#' The soggi function and ChIPprofile object.
#'
#' Manual for soggi and ChIPprofile object
#' 
#' @aliases ChIPprofile ChIPprofile-ChIPprofile soggi
#'
#' @references See \url{http://bioinformatics.csc.mrc.ac.uk} for more details on soGGi workflows
#' @rdname ChIPprofile
#' @docType class
#' @param bamFile Character vector for location of BAM file or bigWig, an rleList or PWM matrix.
#' @param testRanges GRanges object or character vector of BED file location of regions to plot.
#' @param samplename Character vector of sample name. Default is NULL.
#' @param nOfWindows Number of windows to bin regions into for coverage calculations (Default 100)
#' @param FragmentLength Integer vector Predicted or expected fragment length.
#' @param style "Point" for per base pair plot, "percentOfRegion" for normalised length and "region" for combined plot
#' @param distanceAround Distance around centre of region to be used for plotting
#' @param distanceUp Distance upstream from centre of region to be used for plotting
#' @param distanceDown Distance downstream from centre of region to be used for plotting
#' @param distanceInRegionStart Distance into region start 
#' (5' for Watson/positive strand or notspecified strand Regions,3' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceOutRegionStart Distance out from region start 
#' (5' for Watson/positive strand or notspecified strand Regions,3' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceInRegionEnd Distance into region end 
#' (3' for Watson/positive strand or notspecified strand Regions,5' for Crick/negatie strand regions) 
#' for plotting.
#' @param distanceOutRegionEnd Distance out from region end 
#' (3' for Watson/positive strand or notspecified strand Regions,5' for Crick/negatie strand regions) 
#' for plotting.
#' @param paired Is data paired end 
#' @param normalize Calculate coverage as RPM. Presently only RPM available.
#' @param plotBy Score to be used for plotting. Presently only coverage.
#' @param removeDup Remove duplicates before calculating coverage.
#' @param verbose TRUE or FALSE
#' @param format character vector of "BAM", "BigWig", "RleList" or "PWM"
#' @param seqlengths Chromosomes to be used. If missing will report all.
#' @param forceFragment Centre fragment and force consistent fragment width.
#' @param method Character vector of value "bp","bin" or "spline". 
#' The bin method divides a region of interest into equal sized bins of number specified in nOfWindows.
#' Coverage or counts are then summarised within these windows.
#' The spline method creates a spline with the number of spline points as specified in nOfWindows argument.
#' @param downSample Down sample BAM reads to this proportion of orginal.
#' @param genome BSGenome object to be used when using PWM input.
#' @param cutoff Cut-off for idnetifying motifs when using PWM input.
#' @param minFragmentLength Remove fragments smaller than this.
#' @param maxFragmentLength Remove fragments larger than this. 
#' @return ChIPprofile A ChIPprofile object. 
#' @examples
#' data(chipExampleBig)
#' chipExampleBig
#' @export
setClass("ChIPprofile",contains = "RangedSummarizedExperiment",
         slots=c(params="list"
         ))


#' Normalise quantile
#'
#' Quantile normalisation across bins/regions.
#' 
#' @usage
#' \S4method{normaliseQuantiles}{ChIPprofile}(object)
#'
#' @docType methods
#' @name normaliseQuantiles
#' @rdname normaliseQuantiles
#' @aliases normaliseQuantiles normaliseQuantiles,ChIPprofile-method
#' 
#' @author Thomas Carroll
#' @examples
#' data(chipExampleBig)
#' normaliseQuantiles(chipExampleBig)
#' @export
#' @param object A ChIPprofile object
#' @return  A ChIPprofile object containing normalised data
normaliseQuantiles.ChIPprofile <-  function (object)
          {
            
            tempT <- do.call(cbind,lapply(assays(object),function(x)as.vector(x)))
            l <- 1
            
            for(j in 1:ncol(object)){
              tempT[l:(l+nrow(object)-1),] <- normalize.quantiles(tempT[l:(l+nrow(object)-1),])
              l <- l+nrow(object)
            }
            
            qnormAssays <- lapply(1:ncol(tempT),function(x) matrix(tempT[,x],nrow=nrow(object),byrow = FALSE))
            for(c in 1:length(qnormAssays)){
              colnames(qnormAssays[[c]]) <- colnames(assays(object)[[c]])
            }     
            normProfile <- SummarizedExperiment(qnormAssays,rowRanges=rowRanges(object))
            metadata(normProfile) <- metadata(object)
            return(new("ChIPprofile",normProfile,params=object@params))
          }  


setGeneric("normaliseQuantiles", function(object="ChIPprofile") standardGeneric("normaliseQuantiles"))

#' @rdname normaliseQuantiles
#' @export
setMethod("normaliseQuantiles", signature(object="ChIPprofile"), normaliseQuantiles.ChIPprofile)


#' Join, subset and manipulate ChIPprofile objects
#' @rdname manipulateObjects
#' @examples
#' data(chipExampleBig)
#' x <- c(chipExampleBig[[1]],chipExampleBig[[2]])
#' y <- rbind(chipExampleBig[[1]],chipExampleBig[[2]])
#' @return  A ChIPprofile object
#' @export
setMethod("c", "ChIPprofile",
          function (x,...)
          {
            assayList <- lapply(list(x,...),function(x)assays(x)[[1]])
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- unlist(lapply(list(x,...),function(x)metadata(x)$name))
            metadata(subsetProfile)$AlignedReadsInBam <- unlist(lapply(list(x,...),function(x)metadata(x)$AlignedReadsInBam))
            return(new("ChIPprofile",subsetProfile,params=x@params))
            
          }
)

#' @rdname manipulateObjects
#' @export
setMethod("rbind", "ChIPprofile",
          function (x,...,deparse.level=1)
          {
            assayList <- vector("list",length=length(assays(x)))
            regions <- list(x,...)
            for(a in 1:length(assays(x))){
              listTemp <- vector("list",length=length(regions))
              for(r in 1:length(regions)){
                listTemp[[r]] <- assays(regions[[r]])[[a]]
              }
              assayList[[a]] <- do.call(rbind,listTemp)
            }
            newrowRanges <- unlist(GRangesList(lapply(regions,function(x)rowRanges(x))))
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=newrowRanges)
            metadata(subsetProfile)$names <- metadata(x)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x@params))            
          }
)

#' @rdname manipulateObjects
#' @export
setMethod("cbind", "ChIPprofile",
          function (x,...,deparse.level=1)
          {
            assayList <- vector("list",length=length(assays(x)))
            regions <- list(x,...)
            for(a in 1:length(assays(x))){
              listTemp <- vector("list",length=length(regions))
              for(r in 1:length(regions)){
                listTemp[[r]] <- assays(regions[[r]])[[a]]
              }
              assayList[[a]] <- do.call(cbind,listTemp)
            }
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- metadata(x)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x@params))            
          }
)

#' @rdname manipulateObjects
#' @param j Should be missing
#' @export
setMethod("[[", c("ChIPprofile", "ANY", "missing"),
          function(x, i, j, ...)
          {
            subsetProfile <- SummarizedExperiment(assays(x)[[i, ...]],rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- metadata(x)$names[i]
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=x@params))                        
          })


#' @rdname manipulateObjects
#' @export
setMethod("$", "ChIPprofile",
          function(x, name)
          {
            x[[which(metadata(x)$names %in% name)]]
          })

# setMethod("Ops", signature(e1="ChIPprofile", e2="ChIPprofile"),
#           function(e1, e2) {
#             callGeneric(assays(e1)[[1]],assays(e2)[[1]])
#             assaysList <- vector("list",length=length(assays(e1)))
#             for(i in 1:length(assays(e1))){
#               assaysList[[i]] <- callGeneric(assays(e1)[[i]],assays(e2)[[i]]) 
#             }
#             return(assaysList)
#           }
# )

#' Arithmetic operations
#' @rdname Ops
#' @param e1 ChIPprofile object
#' @param e2 ChIPprofile object
#' @examples
#' data(chipExampleBig)
#' chipExampleBig[[1]] + chipExampleBig[[2]]
#' @return  A ChIPprofile object of result of arithmetic operation.
#' @export
setMethod("Ops", signature(e1="ChIPprofile", e2="ChIPprofile"),
          function(e1, e2) {
            assayList <- vector("list",length=length(assays(e1)))
            for(i in 1:length(assays(e1))){
              assayList[[i]] <- callGeneric(assays(e1)[[i]],assays(e2)[[i]]) 
            }
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(e1))
            metadata(subsetProfile)$names <- paste0(metadata(e1)$names,".",metadata(e2)$names)
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(e1)$AlignedReadsInBam[i]
            return(new("ChIPprofile",subsetProfile,params=e1@params))                                    
          }
)

#' @rdname Ops
#' @export
setMethod("Ops", signature(e1="ChIPprofile", e2="numeric"),
          function(e1, e2) {
            assayList <- vector("list",length=length(assays(e1)))
            for(i in 1:length(assays(e1))){
              assayList[[i]] <- callGeneric(assays(e1)[[i]],e2) 
            }
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(e1))
            metadata(subsetProfile)$names <- metadata(e1)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(e1)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=e1@params))                                    
          }
)

#' @rdname Ops
#' @export
setMethod("Ops", signature(e1="numeric", e2="ChIPprofile"),
          function(e1, e2) {
            assayList <- vector("list",length=length(assays(e2)))
            for(i in 1:length(assays(e1))){
              assayList[[i]] <- callGeneric(e1,assays(e2)[[i]]) 
            }
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(e1))
            metadata(subsetProfile)$names <- metadata(e2)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(e2)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=e2@params))                                    
          }
)

#' @rdname Ops
#' @export
setMethod("mean", "ChIPprofile",
          function(x, ...){
            if(missing(...)){
              assayList <- Reduce("+",assays(x))/length(assays(x))
              subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))
              metadata(subsetProfile)$names <- paste0("MeanOf_",metadata(x)$name,collapse="&")
              metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam
              return(new("ChIPprofile",subsetProfile,params=x@params))                                                  
            }else{
              x <- list(x,...)
              assayList <- vector("list",length=length(assays(x[[1]])))
              nameList <- vector("list",length=length(assays(x[[1]])))
              for(a in 1:length(assays(x[[1]]))){
                listTemp <- vector("list",length=length(x))
                nameVec <- vector("character")
                for(r in 1:length(x)){
                  listTemp[[r]] <- assays(x[[r]])[[a]]
                  nameVec <- c(nameVec,metadata(x[[r]])$names[[a]])
                }
                assayList[[a]] <- Reduce("+",listTemp)/length(x)
                nameList[[a]] <- paste0(nameVec,collapse="&")
              }
                subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x[[1]]))
                metadata(subsetProfile)$names <- paste0("MeanOf_",unlist(nameList))
                metadata(subsetProfile)$AlignedReadsInBam <- unlist(lapply(x,function(x)metadata(x)$AlignedReadsInBam))
                return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
                
              }
            }
          
)

#' @rdname Ops
#' @export
setMethod("log2", "ChIPprofile",
          function(x){
            x <- zeroToMin2(x)
            assayList <- lapply(assays(x),log2)
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- metadata(x)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
            
          }
)

#' @rdname Ops
#' @export
setMethod("log", "ChIPprofile",
          function(x,base=exp(1)){
            x <- zeroToMin2(x)
            assayList <- lapply(assays(x),function(x)log(x,base))
            subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(x))
            metadata(subsetProfile)$names <- metadata(x)$names
            metadata(subsetProfile)$AlignedReadsInBam <- metadata(x)$AlignedReadsInBam
            return(new("ChIPprofile",subsetProfile,params=x[[1]]@params))                                                  
            
          }
)


zeroToMin <- function(x){
  for(r in 1:nrow(x)){
    
    temp[r,temp[r,] == 0] <- min(temp[r,temp[r,] != 0])
  }
  return(temp)
}

zeroToMin2 <- function(x){
  for(a in 1:length(assays(x))){
    temp <- assays(x)[[a]]
    ZeroRows <- rowSums(temp) == 0
    temp[ZeroRows,] <- min(temp[temp != 0])
    for(r in 1:nrow(temp)){
      temp[r,temp[r,] == 0] <- min(temp[r,temp[r,] != 0])
    }
    for(r in 1:nrow(temp)){
      temp[r,temp[r,] == 0] <- min(temp[r,temp[r,] != 0])
    }
    assays(x)[[a]] <- temp
  }
  return(x)                                                  
}


#' Normalise ChIPprofiles
#'
#' Various normalisation methods for ChIPprofile objects
#' 
#' @usage
#' \S4method{normalise}{ChIPprofile}(object)
#'
#' @docType methods
#' @name normalise
#' @rdname normalise
#' @aliases normalise normalise,ChIPprofile-method
#' 
#' @author Thomas Carroll
#' @export
#' @param object A ChIPprofile object 
#' @param method A character vector specifying normalisation method. 
#' Currently "rpm" for normalising signal for BAM to total reads, 
#' "quantile" to quantile normalise across samples, 
#' "signalInRegion" to normalise to proportion of signal within intervals,
#' "normaliseSample" to normalise across samples and "normaliseRegions" to apply a normalisation across intervals.
#' @param normFactors A numeric vector used to scale columns or rows.
#' @return  A ChIPprofile object
normalise.ChIPprofile <-  function (object,method="rpm",normFactors=NULL)
{
  
  assaylist <- assays(object)
  if(method=="rpm" & object@params$format == "bam"){
    assayList <- lapply(1:length(assays(object)),
                           function(x) assays(object)[[x]]*metadata(object)$Aligned[x]
                        )
    subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(object))
    metadata(subsetProfile)$names <- metadata(object)$names
    metadata(subsetProfile)$AlignedReadsInBam <- metadata(object)$AlignedReadsInBam
    return(new("ChIPprofile",subsetProfile,params=object[[1]]@params))                                                  
    
    return(subsetProfile)
  }
  if(method=="quantile" ){
    return(normaliseQuantiles(object))
  }
  if(method=="normaliseSamples" & is.numeric(normFactors)){
    return(object*normFactors)
  }
  if(method=="normaliseRegions" & is.numeric(normFactors)){
    assayList <- assays(object)
    for(k in 1:length(assayList)){
      assayList[[k]] <- t(t(assayList[[k]]) * normFactors)
    }
    subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(object))
    metadata(subsetProfile)$names <- metadata(object)$names
    metadata(subsetProfile)$AlignedReadsInBam <- metadata(object)$AlignedReadsInBam
    return(new("ChIPprofile",subsetProfile,params=object[[1]]@params))                                                  
  }
  if(method=="signalInRegion"){
    assayList <- assays(object)
    for(k in 1:length(assayList)){
      assayList[[k]] <- t(t(assayList[[k]])/rowSums(assayList[[k]],na.rm=TRUE))
    }
    subsetProfile <- SummarizedExperiment(assayList,rowRanges=rowRanges(object))
    metadata(subsetProfile)$names <- metadata(object)$names
    metadata(subsetProfile)$AlignedReadsInBam <- metadata(object)$AlignedReadsInBam
    return(new("ChIPprofile",subsetProfile,params=object[[1]]@params))                                                  
  }  
  }
  

setGeneric("normalise", function(object="ChIPprofile",method="rpm",normFactors=NULL) standardGeneric("normalise"))

#' @rdname normalise
#' @examples
#' data(chipExampleBig)
#' normalise(chipExampleBig,method="quantile",normFactors=1)
#' @export
setMethod("normalise", signature(object="ChIPprofile",method="character",normFactors="numeric"), normalise.ChIPprofile)

#' Example ChIPprofiles
#'
#' This dataset contains peaks from ChIP-signal over genes
#'
#' \itemize{
#' \item ChIPprofiles
#' }
#'
#' @docType data
#' @keywords datasets
#' @name chipExampleBig
#' @usage data(chipExampleBig)
#' @return A ChIPprofile object
NULL

#' Example Ikaros signal over peaksets
#'
#' This dataset contains signal over peaks from Ikaros ChIP by two antibodies 
#'
#' \itemize{
#' \item ik_Profiles
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ik_Profiles
#' @usage data(ik_Profiles)
#' @return A ChIPprofile object
NULL

#' Example Ikaros peaksets
#'
#' This dataset contains peaks from Ikaros ChIP by two antibodies 
#'
#' \itemize{
#' \item Ikpeaksets
#' }
#'
#' @docType data
#' @keywords datasets
#' @name ik_Example
#' @usage data(ik_Example)
#' @return A list containing two GRanges objects
NULL

#' Example motif coverage
#'
#' This dataset contains an rlelist of motif coverage
#'
#' \itemize{
#' \item pwmCov
#' }
#'
#' @docType data
#' @keywords datasets
#' @name pwmCov
#' @usage data(pwmCov)
#' @return A rlelist of motif coverage
NULL

#' A single GRange
#'
#' This dataset contains an rlelist of motif coverage
#'
#' \itemize{
#' \item singleGRange
#' }
#'
#' @docType data
#' @keywords datasets
#' @name singleGRange
#' @usage data(singleGRange)
#' @return A single GRanges used in motif coverage example/
NULL
