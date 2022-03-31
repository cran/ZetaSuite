#' Generation of Zeta Plot.
#'
#' A zeta plot is generated from the input Z-score matrix. Zeta plot labels: x-axis: Z-score cutoffs, y-axis: the percentage of readouts that survived at a given Z-score cutoff over the total scored readouts. In order to generate this plot, the range of Z-scores is determined by ranking the absolute value of Zij (Z-score value in row i and column j) from the smallest to the largest. Z-cutoffs next are selected in the range of (-|Znxmx0.9999|, -2) to (2, |Znxmx0.9999|) to excluded the insignificant changes that may result from experimental noise( |Z| < 2, which equals to p-value >0.05). Then, for all Zij within the selected range (both positive range and negative range), the range is divided equally into x bins (the recommended input of x is 100). Thus, the percentage of readouts scored above the Z-cutoff in each bin is determined.
#'
#' @param ZscoreVal zscore value
#'
#' @param negGene negative control dataset, the siRNAs/genes used as negative controls in screening.
#'
#' @param posGene positive control dataset, the siRNAs/genes used as positive controls in screening.
#'
#' @param binNum bin number
#'
#' @param combine combine two direction zeta together(TRUE or FALSE),default FALSE
#'
#' @return A list of data.frames and plots, the data.frame includes 'ZseqList', 'EC_N_I', 'EC_N_D', 'EC_P_I' and 'EC_P_D'. The plot 'EC_jitter_D' and 'EC_jitter_I' are the zeta plot for positive and negative samples.'ZseqList', 'EC_N_I', 'EC_N_D', 'EC_P_I' and 'EC_P_D' are the inputfiles for zeta plot and SVM.R. ZseqList describs the bin size in the zeta plot.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(countMat)
#' data(negGene)
#' data(posGene)
#' ZscoreVal <- Zscore(countMat,negGene)
#' ECList <- EventCoverage(ZscoreVal,negGene,posGene,binNum=100,combine=TRUE)
#'
#' @keywords ZetaSuite
#'
#' @import ggplot2 reshape2 RColorBrewer
#'
#' @importFrom grDevices colorRampPalette dev.off pdf png
#'
#' @importFrom stats quantile
#'
#' @export EventCoverage

EventCoverage <- function (ZscoreVal, negGene, posGene, binNum, combine = TRUE){
    ZscoreVal[is.na(ZscoreVal)] <- 0
	old <- options()
	on.exit(options(old))
    options(digits = 15)
    meltdata <- melt(ZscoreVal)
    if (binNum <= 0){
        stop("binNum should be more than 0")
    }
    if (combine == TRUE){
        minV = quantile(meltdata$value, probs = c(1e-05))
        maxV = quantile(meltdata$value, probs = c(0.99999))
        maxValAbs <- max(abs(minV), abs(maxV))
        stepmax = maxValAbs/binNum
        Zseq_D <- seq(maxValAbs * (-1), -1.3, abs(stepmax))
        Zseq_I <- seq(1.3, maxValAbs, abs(stepmax))
    }else {
        minV = quantile(meltdata$value, probs = c(1e-05))
        stepmin = minV/binNum
        Zseq_D <- seq(minV, 0, abs(stepmin))
        maxV = quantile(meltdata$value, probs = c(0.99999))
        stepmax = maxV/binNum
        Zseq_I <- seq(0, maxV, abs(stepmax))
    }
    ZseqDI <- data.frame(Zseq_D, Zseq_I,stringsAsFactors = FALSE)

    nColCM <- ncol(ZscoreVal)
    EC_D <- matrix(NA, nrow(ZscoreVal), length(Zseq_D))
    colnames(EC_D) <- Zseq_D
    rownames(EC_D) <- rownames(ZscoreVal)
    EC_I <- EC_D
    colnames(EC_I) <- Zseq_I
    for (i in 1:length(Zseq_D)){
        EC_D[, i] <- rowSums(ZscoreVal < Zseq_D[i])/nColCM
        EC_I[, i] <- rowSums(ZscoreVal > Zseq_I[i])/nColCM
    }

    EC_I_N <- EC_I[rownames(EC_I) %in% negGene[, 1], ]
    EC_D_N <- EC_D[rownames(EC_D) %in% negGene[, 1], ]
    EC_I_P <- EC_I[rownames(EC_I) %in% posGene[, 1], ]
    EC_D_P <- EC_D[rownames(EC_D) %in% posGene[, 1], ]
    ECDN <- melt(EC_D_N)
    ECDP <- melt(EC_D_P)
    ECIN <- melt(EC_I_N)
    ECIP <- melt(EC_I_P)

    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    p_EC_jitter_D <- ggplot() + geom_jitter(aes(x = ECDN[, 2], y = ECDN[, 3]), colour = "#67a9cf", size = 0.1) + geom_jitter(aes(x = ECDP[, 2], y = ECDP[, 3]), colour = "#ef8a62", size = 0.1) + ylab("Events percent") + xlab("Zscore") + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values = getPalette(22)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    p_EC_jitter_I <- ggplot() + geom_jitter(aes(x = ECIN[, 2], y = ECIN[, 3]), colour = "#67a9cf", size = 0.1) + geom_jitter(aes(x = ECIP[, 2], y = ECIP[,3]), colour = "#ef8a62", size = 0.1) + ylab("Events percent") + xlab("Zscore") + theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values = getPalette(22)) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ECdataList <- list()
    ECdataList["ZseqList"] <- list(ZseqDI)
    ECdataList["EC_N_I"] <- list(EC_I_N)
    ECdataList["EC_N_D"] <- list(EC_D_N)
    ECdataList["EC_P_I"] <- list(EC_I_P)
    ECdataList["EC_P_D"] <- list(EC_D_P)
    ECplotList <- list()
    ECplotList['EC_jitter_D'] <- list(p_EC_jitter_D)
    ECplotList['EC_jitter_I'] <- list(p_EC_jitter_I)
    return(list(ECdataList,ECplotList))
}

