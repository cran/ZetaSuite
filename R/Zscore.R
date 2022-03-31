#' Z-transformation for input matrix.
#'
#' In this step, the input matrix is transformed to Z-score matrix.
#' 
#' @param countMat input data set. The siRNA/gene x readouts matrix from HTS2 or large-scale RNAi screens.
#'
#' @param negGene negative control dataset, the siRNAs/genes used as negative controls in screening. Z-transfromation according to thses negative control siRNAs/genes for each readout.
#'
#' @details The initial input matrix is arranged in N x M dimension, where each row contains individual functional readouts against a siRNA pool and each column corresponds to individually siRNA pools tested on a given functional readout. Readouts in each column may be thus considered as the data from one-dimensional screen (many-to-one), and thus, the typical Z statistic can be used to evaluate the relative function of individual genes in such column. The conversion is repeated on all columns, thereby converting the raw activity matrix into a matrix. Suppose Nij are the values in the original matrix i (1 <= i <= N siRNA pool) row and j ( 1 <= j <= M readout) column, then Zij = (Nij - uj) / sigma(j), where uj and sigma(j) are the mean and standard deviation of negative control samples in column j.
#' 
#' @return A Z-transformated matrix, where each row represents each knocking-down condition and each column is a specific readout (AS event). The values in the matrix are the normalized values(Z-scores).
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMat)
#' data(negGene)
#' ZscoreVal <- Zscore(countMat,negGene)
#' ZscoreVal[1:5,1:5]
#'
#' @keywords Zscore
#'
#' @importFrom stats sd
#'
#' @export Zscore

Zscore <- function (countMat, negGene) {
    countMatNeg <- countMat[rownames(countMat) %in% negGene[, 1], ]
    countMatNeg$Type <- rep("Negative", nrow(countMatNeg))
    Nmix <- countMatNeg
    rownames(Nmix) <- rownames(countMatNeg)
    matrixData <- countMat
    rownames(matrixData) <- rownames(countMat)
    outputdata <- matrix(0, nrow(matrixData), ncol(Nmix) - 1)
    for (i in 1:(ncol(Nmix) - 1)) {
        outputdata[, i] <- (((matrixData[, i]) - mean(Nmix[, i]))/sd(Nmix[, i]))
    }
    colnames(outputdata) <- colnames(Nmix)[-ncol(Nmix)]
    rownames(outputdata) <- rownames(matrixData)
    return(outputdata)
}
