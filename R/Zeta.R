#' Calculation of zeta and weighted zeta score.
#'
#' This step calculates the Zeta Score based on the two curvecs. Firstly, this step provides another curve above the SVM curve to set a value to represent the regulatory function of gene i. Then, the area between the two curves (the one mentioned above and the SVM curve) is calculated as the Zeta score for this gene. Since the graph of the curves is divided into m bins, then the Zeta score can be calculated as the sum of all the bins' areas that exist between the two curves.
#'
#' @param ZscoreVal input file name.
#'
#' @param ZseqList the list of bins.
#'
#' @param SVMcurve SVM curves for decrease and increase direction.###not always use
#'
#' @param SVM do SVM or not, default is FALSE
#'
#' @return A data.frame where zeta values for all tested knockding-down genes including positive and negative controls. The first column is the direction which knockding-down gene will lead to exon inclusion, whereas the second column is the knock-down genes will lead to exon skipping.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(ZseqList)
#' data(SVMcurve)
#' data(countMat)
#' data(negGene)
#' ZscoreVal <- Zscore(countMat,negGene)
#' zetaData <- Zeta(ZscoreVal,ZseqList,SVM=FALSE)
#'
#' @keywords ZetaSuite zeta
#'
#' @export Zeta

Zeta <- function(ZscoreVal,
                 ZseqList,
                 SVMcurve=NULL,
                 SVM=FALSE){

  Zseq_D <- ZseqList$Zseq_D
  Zseq_I <- ZseqList$Zseq_I
  ZscoreVal[is.na(ZscoreVal)] <- 0
  nColZ <- ncol(ZscoreVal)
  outputdata_D <- rep(0,nrow(ZscoreVal))
  outputdata_I <- rep(0,nrow(ZscoreVal))

  if(SVM==TRUE){
    if(is.null(SVMcurve)==TRUE){
      stop("SVMcurve should not be NULL when SVM is TRUE")
    }
    SVMcurveD <- SVMcurve[,1:2]
    SVMcurveI <- SVMcurve[,3:4]
    for (j in 1:(length(Zseq_D)-1)){
      lengthUse_D <- rowSums(ZscoreVal < Zseq_D[j])
      lengthUse_D_add <- rowSums(ZscoreVal < Zseq_D[j+1])
      conD <- (lengthUse_D/nColZ+lengthUse_D_add/nColZ-SVMcurveD[j,2]-SVMcurveD[j+1,2])*Zseq_D[j+1]*(Zseq_D[j]-Zseq_D[j+1])/2
      conD[conD < 0] <- 0
      outputdata_D<-outputdata_D+conD
      lengthUse_I <- rowSums(ZscoreVal > Zseq_I[j])
      lengthUse_I_add <- rowSums(ZscoreVal > Zseq_I[j+1])
      conI <- (lengthUse_I/nColZ+lengthUse_I_add/nColZ-SVMcurveI[j,2]-SVMcurveI[j+1,2])*Zseq_I[j+1]*(Zseq_I[j+1]-Zseq_I[j])/2
      conI[conI <0] <- 0
      outputdata_I <- outputdata_I + conI
    }
  }else{
    for (j in 1:(length(Zseq_D)-1)){
      lengthUse_D <- rowSums(ZscoreVal < Zseq_D[j])
      lengthUse_D_add <- rowSums(ZscoreVal < Zseq_D[j+1])
      outputdata_D <- outputdata_D + (lengthUse_D/nColZ+lengthUse_D_add/nColZ)*Zseq_D[j+1]*(Zseq_D[j]-Zseq_D[j+1])/2
      lengthUse_I <- rowSums(ZscoreVal > Zseq_I[j])
      lengthUse_I_add <- rowSums(ZscoreVal > Zseq_I[j+1])
      outputdata_I <- outputdata_I + (lengthUse_I/nColZ+lengthUse_I_add/nColZ)*Zseq_I[j+1]*(Zseq_I[j+1]-Zseq_I[j])/2
    }
  }
  outputdata_D <- data.frame(outputdata_D,stringsAsFactors = FALSE)
  colnames(outputdata_D)<-c("Zeta_D")
  outputdata_I <- data.frame(outputdata_I,stringsAsFactors = FALSE)
  colnames(outputdata_I)<-c("Zeta_I")

  zetaDat <- cbind(outputdata_D,outputdata_I)
  return(zetaDat)
}


