#' Plot regions
#'
#' A function to plot regions
#' 
#' @usage
#' \S4method{plotRegion}{ChIPprofile}(object,
#' gts,sampleData,groupData,summariseBy,
#' colourBy,lineBy,groupBy,
#' plotregion,outliers,freeScale)
#'
#'
#' @docType methods
#' @name plotRegion
#' @rdname plotRegion
#' @aliases plotRegion plotRegion,ChIPprofile-method
#' 
#' @author Thomas Carroll
#'
#' @param object A ChIPprofile object 
#' @param gts A list of character vectors or GRangesList
#' @param plotregion region to plot. For combined plots with style "region", may be "start" or "end" to show full resolution of plot of edges.
#' @param groupData Dataframe of metadata for groups
#' @param sampleData Dataframe of metadata for sample
#' @param summariseBy Column names from GRanges elementmetadata. Formula or character vector of column names to use
#' to collapse genomic ranges to summarised profiles.
#' summariseBy can not be used injustion with groups specified by gts argument.
#' @param colourBy Character vector or formula of either column names from colData(object) containing
#' sample metadata or character vector "group" to colour by groups in gts
#' @param lineBy Character vector or formula of either column names from colData(object) containing
#' sample metadata or character vector "group" to set line type by groups in gts
#' @param groupBy Character vector or formula of either column names from colData(object) containing
#' sample metadata or character "group" to colour by groups in gts
#' @param outliers A numeric vector of length 1 containing proportion from limits to windsorise.]
#' @param freeScale TRUE or FALSE to set whether ggplot 2 facets have their own scales. 
#' Useful for comparing multiple samples of differing depths without normalisation. Default is FALSE.
#' @return  A gg object from ggplot2
#' @examples
#' data(chipExampleBig)
#' plotRegion(chipExampleBig[[2]])
plotRegion.ChIPprofile <- function(object,gts=NULL,sampleData=NULL,groupData=NULL,summariseBy=NULL,colourBy=NULL,lineBy=NULL,groupBy=NULL,plotregion="full",outliers=NULL,freeScale=FALSE)
{
  
## This function can work using two main grouping logics
## When gts is supplied, groupBy, colourBy, lineBy work only with columns
## from the sample metadata found in colData(object).
## If the gts argument is specified then summariseBy must be single character vector
##  specifies which metadata column which column to be used for grouping
##
## When summariseBy is specified as formula or character vector longer  
##  than one then maintain as/coerced into a formula and gts is ignored.
## groupBy, colourBy, lineBy must be parts of formula as assessed by terms()  
## 
## gts may be named list of character vectors referring to rownames(object) or a granges
## When 
  #app <- lapply(gsets,function(x){colMeans(assays(object)[[1]][rowRanges(object)$name %in% x,])})
  nOfWindows <- object@params$nOfWindows
  
  ## When running with gts option
  ## SummariseBy can now only be used to select column for gts to be match to
  ## and colourBy, lineBy and groupBy will refer to sample metadata or be "group"
  
  if(!is.null(gts)){
    
    profileList <- list()
    
    ## Start cycling through assays
    for(p in 1:length(assays(object))){
    
    ## extract profile matrix
      profileTemp <- assays(object)[[p]]
    
    ## profiles are subset using subsetProfile.
    ## this allows for subsets by lists of character vector 
    ## or granges lists.colmeans is then run on these subsets
    ## Alternatively windsoring (see method) of subsets then colmeans.
      
      if(!is.null(outliers)){
        profileTempList <- lapply(gts,function(x)
          colMeans(winsorizeMatrix(subsetProfile(profileTemp,x,rowRanges(object),summariseBy),outliers,1-outliers))
          )         
      }else{
        profileTempList <- lapply(gts,function(x)colMeans(subsetProfile(profileTemp,x,rowRanges(object),summariseBy))) 
      }
    
    ## Create melted data frame for ggplot and attach index
    
      profileMatTemp <- melt(as.data.frame(do.call(cbind,profileTempList)))
    
    ## and attach index for different styles of plots
    if(object@params$style=="region" & plotregion=="full"){
        axisIndex=c(seq(1,(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)),
                    (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                    (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                      seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)))
      }
      if(object@params$style=="point"){
        axisIndex=c(seq(1,(object@params$distanceUp+object@params$distanceDown+1)))
      }
      if(object@params$style=="percentOfRegion"){
        axisIndex=c(seq(1,((nOfWindows*((object@params$distanceAround)/100))*2)+nOfWindows))
      } 
    
    ## Add Sample name, group name and index to dataframe
      profileFrame <-data.frame("xIndex"=axisIndex,Group=profileMatTemp[,1],Sample=basename(unlist(metadata(object)["names"]))[p],Score=profileMatTemp[,2])
      
      profileList[[p]] <- profileFrame
    }  
    
    ## Join profiles from each sample and fix column names
    meltedProfileFrame <- do.call(rbind,profileList)
    colnames(meltedProfileFrame) <- c("xIndex","Group","Sample","Score")

  }else if(is.null(gts) & !is.null(summariseBy)){
      profileList <- list()
      
      ## Start cycling through assays
      for(p in 1:length(assays(object))){
            
            ## extract profile matrix
          profileTemp <- assays(object)[[p]]
            
          mcols(rowRanges(object))$summariseCol <- apply(as.data.frame(mcols(rowRanges(object))[,summariseBy,drop=FALSE]) , 1 , paste , collapse = "-" )
          gts <- as.list(unique(mcols(rowRanges(object))$summariseCol))
          summariseBy <- "summariseCol"
          names(gts) <- unlist(gts)
          if(!is.null(outliers)){
            profileTempList <- lapply(gts,function(x)
              colMeans(winsorizeMatrix(subsetProfile(profileTemp,x,rowRanges(object),summariseBy),outliers,1-outliers))
            )         
          }else{
            profileTempList <- lapply(gts,function(x)colMeans(subsetProfile(profileTemp,x,rowRanges(object),summariseBy))) 
          }
      
    
    ## Create melted data frame for ggplot and attach index
    
    profileMatTemp <- melt(as.data.frame(do.call(cbind,profileTempList)))
    if(object@params$style=="region" & plotregion=="full"){
      axisIndex=c(seq(1,(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)),
                  (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                  (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                    seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)))
    }
    if(object@params$style=="point"){
      axisIndex=c(seq(1,(object@params$distanceUp+object@params$distanceDown+1)))
    }
    if(object@params$style=="percentOfRegion"){
      axisIndex=c(seq(1,((nOfWindows*((object@params$distanceAround)/100))*2)+nOfWindows))
    } 
    
    ## Add Sample name, group name and index to dataframe
    profileFrame <-data.frame("xIndex"=axisIndex,Group=profileMatTemp[,1],Sample=basename(unlist(metadata(object)["names"]))[p],Score=profileMatTemp[,2])
    
    profileList[[p]] <- profileFrame
  }  
  
  ## Join profiles from each sample and fix column names
  meltedProfileFrame <- do.call(rbind,profileList)
  colnames(meltedProfileFrame) <- c("xIndex","Group","Sample","Score")
  
  }else{
    
    ## When no summariseBy or gts supplied colmeans or
    ## windsorised colmeans of whole profile matrix
    
    if(!is.null(outliers)){
      profileList <- lapply(c(assays(object)),function(x)
        colMeans(winsorizeMatrix(x,outliers,1-outliers))
      )         
    }else{    
      profileList <- lapply(c(assays(object)),colMeans)
    }
    
    ## Join multiple assays/samples
    
    profileFrame <- do.call(cbind,profileList)
    colnames(profileFrame) <- basename(unlist(metadata(object)["names"]))

    ## Attach index for different styles of plots
    
    if(object@params$style=="region"){    
      axisIndex=c(seq(1,(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)),
                  (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+seq(1,object@params$nOfWindows)*100,
                  (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+
                    seq(1,(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)))
    }
    if(object@params$style=="point"){
      axisIndex=c(seq(1,(object@params$distanceUp+object@params$distanceDown+1)))
    }
    if(object@params$style=="percentOfRegion"){
      axisIndex=c(seq(1,((object@params$nOfWindows*((object@params$distanceAround)/100))*2)+object@params$nOfWindows))
    }     
    rownames(profileFrame) <- axisIndex
    meltedProfileFrame <- melt(profileFrame)
    colnames(meltedProfileFrame) <- c("xIndex","Sample","Score")
  }
  #profileList <- lapply(c(assays(object)),function(y)lapply(gsets,function(x){colMeans(y[rowRanges(object)$name %in% x,])}))
  
  ## Create geom_path plot
  if(!is.null(gts) & !is.null(groupData)){
    meltedProfileFrame <- merge(meltedProfileFrame,groupData,by.x=2,by.y=1,all.x=TRUE,all.y=FALSE,sort=FALSE)
  }
  if(!is.null(sampleData)){
    sampleNameCol <- which(colnames(meltedProfileFrame) %in% "Sample")
    meltedProfileFrame <- merge(meltedProfileFrame,sampleData,by.x=sampleNameCol,by.y=1,all.x=TRUE,all.y=FALSE,sort=FALSE)
  }
  
  P <- ggplot(meltedProfileFrame,
              aes_string(x="xIndex",y="Score"))+geom_path(alpha = 1,size=1.3)+xlim(0,max(axisIndex))+ylab("Score")+theme(axis.title.y=element_text(angle=0))
  
  ## Add scales depending on style and region being plotted
  
  if(object@params$style=="region" & plotregion =="full"){
    P <- P + scale_x_continuous(breaks=c(1,object@params$distanceOutRegionStart+1,
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(25/100)*(object@params$nOfWindows*100),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(75/100)*(object@params$nOfWindows*100),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100),
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+object@params$distanceInRegionEnd+1,
                                       (object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+(object@params$distanceInRegionEnd+object@params$distanceOutRegionEnd+1)),
                              labels=c(paste0("Start-",object@params$distanceOutRegionStart),
                                       "Start",
                                       paste0("Start+",object@params$distanceInRegionStart),
                                       "25%",
                                       "75%",
                                       paste0("End-",object@params$distanceInRegionEnd),
                                       "End",
                                       paste0("End+",object@params$distanceOutRegionEnd)
                                       )
                                )+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }
  if(object@params$style=="point"){
    P <- P + scale_x_continuous(breaks=c(1,object@params$distanceUp+1,object@params$distanceUp+1+object@params$distanceDown),
                              labels=c(paste0("Centre-",object@params$distanceUp),"Centre",paste0("Centre+",object@params$distanceDown)))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }
  if(object@params$style=="percentOfRegion"){
    P <- P + scale_x_continuous(breaks=c(1,(nOfWindows*((object@params$distanceAround)/100)),
                                            (nOfWindows*((object@params$distanceAround)/100))+nOfWindows,
                                             2*(nOfWindows*((object@params$distanceAround)/100))+nOfWindows),
                                labels=c(paste0("Start-",(nOfWindows*((object@params$distanceAround)/100)),"%"),
                                                "Start",
                                                "End",
                                                paste0("End+",(nOfWindows*((object@params$distanceAround)/100)),"%")
                                                ))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))
  }  
  if(object@params$style=="region" & plotregion =="start"){
    P <- P + scale_x_continuous(breaks=c(1,object@params$distanceOutRegionStart+1,object@params$distanceInRegionStart+1+object@params$distanceOutRegionStart),
                              labels=c(paste0("Start-",object@params$distanceOutRegionStart),"Centre",paste0("Start-",object@params$distanceInRegionStart)),
                              limits=c(1,object@params$distanceInRegionStart+1+object@params$distanceOutRegionStart))+
      theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))    
  }
  if(object@params$style=="region" & plotregion =="end"){
    P <- P + scale_x_continuous(breaks=(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+c(1,object@params$distanceInRegionEnd+1,object@params$distanceInRegionEnd+1+object@params$distanceOutRegionEnd),
                              labels=c(paste0("Start-",object@params$distanceInRegionEnd),"Centre",paste0("Start-",object@params$distanceOutRegionEnd)),
                              limits=(object@params$distanceOutRegionStart+object@params$distanceInRegionStart+1)+(object@params$nOfWindows*100)+c(1,object@params$distanceInRegionEnd+1+object@params$distanceOutRegionEnd))+
                                theme(axis.text.x  = element_text(angle=45, vjust=0.5, size=12))    
  } 
  
  

  if(!is.null(gts)){
    if(is.null(groupBy) & is.null(colourBy) & is.null(lineBy)){
      groupBy <- "Group"
    }
    #P <- P+aes(group="group")+aes_string(colour=colourBy,linetype=lineBy)   
  }
  if(is.null(groupBy) & is.null(colourBy) & is.null(lineBy)){
    groupBy <- "Sample"
  }

  
  if(!is.null(groupBy)){
    
    facet <- facet_wrap(
      formula(paste("~",paste(groupBy,collapse="+")))
    )
    P <- P + facet
    
  }
  
  if(freeScale){
    P$facet$free$y <- TRUE 
  }
  P <- P+aes_string(colour=colourBy,linetype=lineBy)
  

  
  return(P)
}

winsorizeMatrix <- function(mat,limitlow,limithigh){
  apply(mat,2,function(x)winsorizeVector(x,limitlow,limithigh))
}
winsorizeVector <- function(vect,limitlow,limithigh){
  qs <- quantile(vect,c(limitlow,limithigh))
  vect[vect < qs[1]] <- qs[1]  
  vect[vect > qs[2]] <- qs[2]  
  vect
}
setGeneric("plotRegion", function(object="ChIPprofile",gts=NULL,sampleData=NULL,groupData=NULL,summariseBy=NULL,colourBy=NULL,lineBy=NULL,groupBy=NULL,plotregion="character",outliers=NULL,freeScale=FALSE) standardGeneric("plotRegion"))

#' @rdname plotRegion
#' @export
setMethod("plotRegion", signature(object="ChIPprofile"), plotRegion.ChIPprofile)

subsetProfile <- function(profile,group,granges,summariseColumn){
  if(class(group) == "character"){
    if(is.null(summariseColumn)){
      return(profile[rownames(profile) %in% group,])
    }else{
      return(profile[mcols(granges)[,summariseColumn] %in% group,])
    }
  }
  if(class(group) == "GRanges"){         
    return(profile[granges %over% group,])
  }
}
