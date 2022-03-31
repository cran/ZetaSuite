## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

## ---- eval=F, echo=T----------------------------------------------------------
#  ## install ZetaSuite R package
#  install.packages("ZetaSuite")
#  
#  ## load ZetaSuite R pakcage
#  library(ZetaSuite)

## ----workflow1, message=FALSE, include=FALSE----------------------------------
library(ZetaSuite)
data(countMat)
data(negGene)
data(posGene)
qcList <- QC(countMat,negGene,posGene)

## ----workflow1.1, eval=F, echo=T----------------------------------------------
#  library(ZetaSuite)
#  data(countMat)
#  data(negGene)
#  data(posGene)
#  qcList <- QC(countMat,negGene,posGene)

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
qcList[[2]]

## ---- fig.height = 3, fig.width = 8, fig.align = "center"---------------------
qcList[[1]]

## ---- fig.height = 4, fig.width = 8, fig.align = "center"---------------------
grid::grid.draw(qcList[[3]])

## ---- fig.height = 4, fig.width = 6, fig.align = "center"---------------------
qcList[[4]]

## ----workflow2----------------------------------------------------------------
data(countMat)
data(negGene)
ZscoreVal <- Zscore(countMat,negGene)
ZscoreVal[1:5,1:5]

## ----workflow3, fig.height = 4, fig.width = 4, fig.align = "center"-----------
data(countMat)
data(negGene)
data(posGene)
ZscoreVal <- Zscore(countMat,negGene)
ECList <- EventCoverage(ZscoreVal,negGene,posGene,binNum=100,combine=TRUE)
ECList[[2]][[1]]
ECList[[2]][[2]]

## ----workflow4----------------------------------------------------------------
data(ZseqList)
data(SVMcurve)
data(countMat)
data(negGene)
ZscoreVal <- Zscore(countMat,negGene)
zetaData <- Zeta(ZscoreVal,ZseqList,SVM=FALSE)

## ----workflow5, fig.height = 3, fig.width = 6, fig.align = "center"-----------
data(nonExpGene)
data(negGene)
data(posGene)
data(ZseqList)
data(countMat)
ZscoreVal <- Zscore(countMat,negGene)
zetaData <- Zeta(ZscoreVal,ZseqList,SVM=FALSE)
cutoffval <- FDRcutoff(zetaData,negGene,posGene,nonExpGene,combine=TRUE)
cutoffval[[2]][[1]]
cutoffval[[2]][[2]]

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

