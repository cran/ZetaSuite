#' Find a cut-off according to screen strength.
#'
#' Find a cutoff according to the Screen Strength (SS) and graph the Screen Strength plot.
#' Zeta score is used to rank genes, and then, SS is calculated to define a suitable cutoff so that the cutoff can define hits at different confidence intervals.
#' Formula of SS: SS = 1 - aFDR/bFDR, where aFDR (apparent FDR)  = number of non-expressors identified at hits divided by the total number of hits, bFDR (baseline FDR) = total number of non-expressors divided by all screened genes.
#' SS plot labels: x-axis: zeta score, y-axis: Screen Strength,
#' SS value is determined at each bin (m bin in total), then connect individual SS value to generate a simulated SS curve based on balance points. Users may choose one or multiple balance point as the different SS intervals.
#'
#' @param zetaData ZetaScore file calculated by ZetaSuite.
#'
#' @param negGene negative control dataset, the siRNAs/genes used as negative controls in screening.
#'
#' @param posGene positive control dataset, the siRNAs/genes used as positive controls in screening.
#'
#' @param nonExpGene non-expressed gene
#'
#' @param combine combine two direction zeta together(TRUE or FALSE),default FALSE
#'
#' @return A list of data.frame and plots, the data.frame is cut off matrix with 6 columns including "Cut_Off","aFDR", "SS","TotalHits","Num_nonExp" and "Type". Plots includes 'Zeta_type' and 'SS_cutOff'.
#'
#' @author Yajing Hao, Shuyang Zhang, Junhui Li, Guofeng Zhao, Xiang-Dong Fu
#'
#' @examples
#' data(nonExpGene)
#' data(negGene)
#' data(posGene)
#' data(ZseqList)
#' data(countMat)
#' ZscoreVal <- Zscore(countMat,negGene)
#' zetaData <- Zeta(ZscoreVal,ZseqList,SVM=FALSE)
#' cutoffval <- FDRcutoff(zetaData,negGene,posGene,nonExpGene,combine=TRUE)
#'
#' @keywords ZetaSuite FDR cutoff
#'
#' @import ggplot2
#'
#' @importFrom grDevices dev.off pdf
#'
#' @export FDRcutoff

FDRcutoff <- function (zetaData, negGene, posGene, nonExpGene, combine = FALSE) {
  zetaData$type <- rep("Gene", nrow(zetaData))
  zetaData[rownames(zetaData) %in% negGene[, 1], "type"] <- "NS_mix"
  zetaData[rownames(zetaData) %in% posGene[, 1], "type"] <- "Positive"
  zetaData[rownames(zetaData) %in% nonExpGene[, 1], "type"] <- "non_exp"
  nonGene_D <- zetaData[zetaData$type %in% c("non_exp", "Gene"), ]

  if (combine == FALSE) {
    maxD <- sort(nonGene_D[, "Zeta_D"], decreasing = TRUE)[10]
    maxI <- sort(nonGene_D[, "Zeta_I"], decreasing = TRUE)[10]
    minD <- sort(nonGene_D[, "Zeta_D"], decreasing = FALSE)[1]
    minI <- sort(nonGene_D[, "Zeta_I"], decreasing = FALSE)[1]
    stepD <- (maxD - minD)/100
    stepI <- (maxI - minI)/100
    sum1 <- nrow(nonGene_D)
    sum2 <- sum(nonGene_D[, "type"] %in% "non_exp")
    iFDR <- sum2/(sum1 - 1)
    seqD <- seq(minD, maxD, stepD)
    FDR_cutOff_de <- matrix(NA, length(seqD), 6)
    index <- 1
    for (num in seqD) {
      totalNum <- sum(zetaData[, "Zeta_D"] >= num & zetaData[, "type"] %in% c("non_exp", "Gene"))
      numNexp <- sum(zetaData[, "Zeta_D"] >= num & zetaData[, "type"] %in% c("non_exp"))
      FDR_Nexp <- numNexp/totalNum
      screen_Stress <- (iFDR - FDR_Nexp)/iFDR
      FDR_cutOff_de[index, ] <- c(num, FDR_Nexp,screen_Stress,totalNum, numNexp, "Decrease")
      index <- index + 1
    }
    seqI <- seq(minI, maxI, stepI)
    FDR_cutOff_in <- matrix(NA, length(seqI), 6)
    index <- 1
    for (num in seqI) {
      totalNum <- sum(zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp", "Gene"))
      numNexp <- sum(zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp"))
      FDR_Nexp <- numNexp/totalNum
      screen_Stress <- (iFDR - FDR_Nexp)/iFDR
      FDR_cutOff_in[index, ] <- c(num, FDR_Nexp, screen_Stress, totalNum, numNexp, "Increase")
      index <- index + 1
    }
    FDR_cutOff <- rbind(FDR_cutOff_de, FDR_cutOff_in)
  } else {
    maxD <- sort(nonGene_D[, 1] + nonGene_D[, 2], decreasing = TRUE)[10]
    minD <- sort(nonGene_D[, 1] + nonGene_D[, 2], decreasing = FALSE)[1]
    stepD <- (maxD - minD)/100
    sum1 <- nrow(nonGene_D)
    sum2 <- sum(nonGene_D[, 3] %in% "non_exp")
    iFDR <- sum2/(sum1 - 1)
    seqD <- seq(minD, maxD, stepD)
    FDR_cutOff <- matrix(NA, length(seqD), 6)
    index <- 1
    for (num in seqD) {
      totalNum <- sum(zetaData[, "Zeta_D"] + zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp", "Gene"))
      numNexp <- sum(zetaData[, "Zeta_D"] + zetaData[, "Zeta_I"] >= num & zetaData[, "type"] %in% c("non_exp"))
      FDR_Nexp <- numNexp/totalNum
      screen_Stress <- (iFDR - FDR_Nexp)/iFDR
      FDR_cutOff[index, ] <- c(num, FDR_Nexp, screen_Stress, totalNum, numNexp, "Combine")
      index <- index + 1
    }
  }
  FDR_cutOff <- as.data.frame(FDR_cutOff,stringsAsFactors = FALSE)
  colnames(FDR_cutOff) <- c("Cut_Off","aFDR","SS","TotalHits","Num_nonExp","Type")
  FDR_cutOff$SS <- as.numeric(FDR_cutOff$SS)
  FDR_cutOff$Cut_Off <- as.numeric(FDR_cutOff$Cut_Off)
  zetaData_NS <- zetaData[zetaData$type != "NS_mix", ]

  p1 <- ggplot(zetaData_NS) + geom_jitter(aes_string(x = "type", y = "Zeta_D", col = "type")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + scale_color_manual(values = c("#c994c7", "#67a9cf", "#ef8a62"))
  p2 <- ggplot(zetaData_NS) + geom_jitter(aes_string(x = "type", y = "Zeta_I", col = "type")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + scale_color_manual(values = c("#c994c7", "#67a9cf", "#ef8a62"))
  p_Zeta_type <- gridExtra::grid.arrange(p1, p2, nrow = 1)

  if (combine == FALSE) {
    fdtss <- FDR_cutOff[FDR_cutOff$SS < 0.9, ]
    Dec <- fdtss[fdtss$Type == "Decrease", ]
    Inc <- fdtss[fdtss$Type == "Increase", ]

    p1 <- ggplot(Dec, aes_string(x = "Cut_Off", y = "SS", col = "Type")) + geom_point() + geom_smooth(span = 0.2) + theme_bw() + xlab("Zeta Score") + ylab("Screen strength") + theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    p2 <- ggplot(Inc, aes_string(x = "Cut_Off", y = "SS", col = "Type")) + geom_point() + geom_smooth(span = 0.2) + theme_bw() + xlab("Zeta Score") + ylab("Screen strength") + theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
    p_SS_cutOff <- gridExtra::grid.arrange(p1, p2, nrow = 1)

  } else {
    fdtss <- FDR_cutOff[FDR_cutOff$SS < 1, ]

    p_SS_cutOff <- ggplot(fdtss, aes_string(x = "Cut_Off", y = "SS", col = "Type")) + geom_point() + geom_smooth(span = 0.2) + theme_bw() + xlab("Zeta Score") + ylab("Screen strength") + theme(legend.position = c(0.8, 0.2), legend.title = element_blank())
  }
  plotList <- list()
  plotList['Zeta_type'] <- list(p_Zeta_type)
  plotList['SS_cutOff'] <- list(p_SS_cutOff)
  return(list(FDR_cutOff,plotList))
}
