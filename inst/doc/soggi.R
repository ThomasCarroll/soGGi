## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,dev="png",fig.show="hide",
               fig.width=4,fig.height=4.5,
               message=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----loadSoggi, echo=FALSE-------------------------------------------------
library("soGGi")

## ----include=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----include=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----options, results="hide", echo=FALSE--------------------------------------
options(digits=3, width=80, prompt=" ", continue=" ")

## ----quick1, eval=FALSE , fig.width=10, fig.height=6--------------------------
#  library(soGGi)
#  chipExample <- regionPlot("pathToBAM/mybam.bam",myGRangesObject,format="bam")

## ----quick1_1, eval=TRUE , fig.width=10, fig.height=6-------------------------
library(soGGi)
data(chipExampleBig)
chipExampleBig

## ----quick1_2, eval=TRUE , fig.width=10, fig.height=6-------------------------
chipExampleBig[[1]]
chipExampleBig$highdnase

## ----quick1_3, eval=TRUE , fig.width=10, fig.height=6-------------------------
c(chipExampleBig[[1]],chipExampleBig[[2]])
rbind(chipExampleBig[[1]],chipExampleBig[[2]])

## ----quicka, eval=TRUE , fig.width=10, fig.height=6---------------------------
plotRegion(chipExampleBig[[3]])

## ----quick2, eval=TRUE , fig.width=10, fig.height=6---------------------------
library(ggplot2)
plotRegion(chipExampleBig,colourBy="Sample", groupBy="Sample", freeScale=TRUE)

## ----quick2b, eval=TRUE , fig.width=10, fig.height=6--------------------------
library(ggplot2)
plotRegion(chipExampleBig,colourBy="Sample", outliers=0.01, groupBy="Sample",freeScale=TRUE)

## ----quick3, eval=TRUE , fig.width=10, fig.height=6---------------------------
library(GenomicRanges)
subsetsCharacter <- list(first25 = (as.vector(rowRanges(chipExampleBig[[1]])$name[1:25])), fourth25 = as.vector(rowRanges(chipExampleBig[[1]])$name[76:100]))

subsetsGRanges <- GRangesList(low=(rowRanges(chipExampleBig[[1]])[1:25]), high=rowRanges(chipExampleBig[[2]])[76:100])

plotRegion(chipExampleBig[[1]],gts=subsetsCharacter,summariseBy = "name")
plotRegion(chipExampleBig[[1]],gts=subsetsGRanges)


## ----quick4, eval=TRUE , fig.width=10, fig.height=6---------------------------
pol_Profiles <- c((chipExampleBig$highPol+chipExampleBig$midPol)
, (chipExampleBig$highPol_Rep2+chipExampleBig$midPol_Rep2))
plotRegion(pol_Profiles,colourBy="Sample",outliers=0.01, groupBy="Sample", freeScale=TRUE)

## ----quick6, eval=TRUE , fig.width=10, fig.height=6---------------------------
log2Profiles <- log2(chipExampleBig)
                     
plotRegion(log2Profiles,colourBy="Sample",outliers=0.01, groupBy="Sample",freeScale=TRUE)


## ----quick6b, eval=TRUE , fig.width=10, fig.height=6--------------------------
log2Polhigh <- mean(log2Profiles$highPol, log2Profiles$highPol_Rep2)
log2Polmid <- mean(log2Profiles$midPol, log2Profiles$midPol_Rep2)                     
diffPol <- log2Polhigh-log2Polmid

                     
diffh3k9ac <- log2Profiles$highk9ac-log2Profiles$midk9ac


plotRegion(c(diffPol,diffh3k9ac),colourBy="Sample",outliers=0.01, groupBy="Sample", freeScale=TRUE)


## ----quick7, eval=TRUE , fig.width=10, fig.height=6---------------------------
normHighPol <- normalise(c(chipExampleBig$highPol, chipExampleBig$highPol_Rep2), method="quantile",normFactors = 1)
normMidPol <- normalise(c(chipExampleBig$midPol, chipExampleBig$midPol_Rep2), method="quantile",normFactors = 1)
         
normPol <-c(normHighPol$highPol, normHighPol$highPol_Rep2, normMidPol$midPol, normMidPol$midPol_Rep2)
plotRegion(normPol,colourBy="Sample",outliers=0.01, groupBy="Sample", freeScale=TRUE)


## ----quick8, eval=TRUE , fig.width=10, fig.height=5---------------------------
data(ik_Example)
ik_Example
peakSetCombinations <- groupByOverlaps(ik_Example)
peakSetCombinations

## ----quick10, eval=TRUE , fig.width=10, fig.height=5--------------------------
data(ik_Profiles)
ik_Profiles
log2Ik_Profiles <- log2(ik_Profiles)
plotRegion(log2Ik_Profiles,outliers=0.01,gts=peakSetCombinations, groupBy="Group",colourBy="Sample", freeScale=TRUE)
plotRegion(log2Ik_Profiles[[1]] - log2Ik_Profiles[[2]] ,outliers=0.01,gts=peakSetCombinations, groupBy="Group", colourBy="Sample", freeScale=FALSE)


## ----sessionInfo, results='asis', eval=TRUE-----------------------------------
toLatex(sessionInfo())

