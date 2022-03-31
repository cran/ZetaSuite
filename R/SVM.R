#' Find a SVM curve to separate positive and negative controls.
#'
#' Radical kernel SVM is constructed to maximally separate positive controls from negative controls in the prior defined Z range using e1071 packages of R, and therefore, the SVM curve is generated.
#'
#' @param ECdataList data list of output EventCoverage, names of list shoule be 'EC_N_D', 'EC_P_D', 'EC_N_I', 'EC_P_I' and 'ZseqList'
#'
#' @return A list of data.frame, including 'cutOffD' and 'cutOffI'.cutOffD and cutOffI are the deduced SVM. 
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMat)
#' data(negGene)
#' data(posGene)
#' ZscoreVal <- Zscore(countMat,negGene)
#' ECdataList <- EventCoverage(ZscoreVal,negGene,posGene,binNum=10,combine=TRUE)
#' \donttest{SVM(ECdataList)}
#'
#' @keywords support vector machine
#'
#' @import e1071 ggplot2 reshape2
#'
#' @importFrom stats predict
#'
#' @export SVM

SVM <- function(ECdataList){
  old <- options()
  on.exit(options(old))
  options(digits=15)
  EC_N_D <- melt(ECdataList[[1]][["EC_N_D"]])
  EC_P_D <- melt(ECdataList[[1]][["EC_P_D"]])
  P_D <- cbind(EC_P_D,rep("Positive",nrow(EC_P_D)))
  N_D <- cbind(EC_N_D,rep("Negative",nrow(EC_N_D)))
  P_D <- P_D[,-1]
  N_D <- N_D[,-1]
  colnames(P_D) <- c("Zscore","Percent","label")
  colnames(N_D) <- c("Zscore","Percent","label")
  SVMInput_D <- rbind(P_D,N_D)
  SVMInput_D$label <- as.factor(SVMInput_D$label)

  EC_N_I <- melt(ECdataList[[1]][["EC_N_I"]])
  EC_P_I <- melt(ECdataList[[1]][["EC_P_I"]])
  P_I <- cbind(EC_P_I,rep("Positive",nrow(EC_P_I)))
  N_I <- cbind(EC_N_I,rep("Negative",nrow(EC_N_I)))
  P_I <- P_I[,-1]
  N_I <- N_I[,-1]
  colnames(P_I)<-c("Zscore","Percent","label")
  colnames(N_I)<-c("Zscore","Percent","label")
  SVMInput_I <- rbind(P_I,N_I)
  SVMInput_I$label <- as.factor(SVMInput_I$label)

  svm_D_tune <- svm(label ~ ., data=SVMInput_D, type="C-classification", kernel="radial", cost=20, gamma=3, scale=FALSE,probability=TRUE)
  #svm_model_after_tune
  grid <- expand.grid(seq(min(SVMInput_D[,1]), max(SVMInput_D[,1]),length.out=length(unique(SVMInput_D[,1]))), seq(min(SVMInput_D[,2]), max(SVMInput_D[,2]),length.out=length(unique(SVMInput_D[,1]))))
  names(grid) <- names(SVMInput_D)[1:2]
  preds <- predict(svm_D_tune, grid)
  svm_line_D <- data.frame(grid, preds,stringsAsFactors = FALSE)

  svm_I_tune <- svm(label ~ ., data=SVMInput_I, type="C-classification", kernel="radial", cost=50, gmma=2, scale=FALSE,probability=TRUE)
  grid <- expand.grid(seq(min(SVMInput_I[,1]), max(SVMInput_I[,1]),length.out=length(unique(SVMInput_I[,1]))), seq(min(SVMInput_I[,2]), max(SVMInput_I[,2]),length.out=length(unique(SVMInput_I[,1]))))
  names(grid) <- names(SVMInput_I)[1:2]
  preds <- predict(svm_I_tune, grid)
  svm_line_I <- data.frame(grid, preds,stringsAsFactors = FALSE)

  ZseqList <- ECdataList[[1]][["ZseqList"]]
  lineCutOffD <- NULL
  lineCutOffI <- NULL
  for(j in 1:nrow(ZseqList)){
    numberD <- svm_line_D[,1]-ZseqList[j,1]
    numberI <- svm_line_I[,1]-ZseqList[j,2]
    indexD <- which(numberD <= 0.00000001 & numberD >= -0.00000001)
    indexI <- which(numberI <= 0.00000001 & numberI >= -0.00000001)
    if(length(indexD)>0){
      indexD_neg <- indexD[svm_line_D[indexD,3]=="Negative"]
      indexD_neg_max <- indexD_neg[which.max(svm_line_D[indexD_neg,2])]
      lineCutOffD <- rbind(lineCutOffD,svm_line_D[indexD_neg_max,c(1:2)])
    }
    if(length(indexI)>0){
      indexI_neg <- indexI[svm_line_I[indexI,3]=="Negative"]
      indexI_neg_max <- indexI_neg[which.max(svm_line_I[indexI_neg,2])]
      lineCutOffI <- rbind(lineCutOffI,svm_line_I[indexI_neg_max,c(1:2)])
    }
  }
  svmList <- list()
  svmList["cutOffD"] <- list(lineCutOffD)
  svmList["cutOffI"] <- list(lineCutOffI)
  return(svmList)
}
