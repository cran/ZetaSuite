\name{countMat}
\alias{countMat}
\docType{data}
\title{
  Subsampled data from in-house HTS2 screening for global splicing regulators.
}
\description{
  A data frame with 1609 individual screened genes and 100 functional readouts. The data was generated from a siRNA screen for global splicing regulators. In this screen, we interrogated ~400 endogenous alternative splicing (AS) events by using an oligo ligation-based strategy to quantify 18,480 pools of siRNAs against annotated protein-coding genes in the human genome.
}
\usage{data("countMat")}
\format{
  A data frame with 1609 observations on the following 100 variables
  \describe{
    A data frame with 1609 observations on the following 100 maker variables.Each row represents gene with specific knocking-down siRNA pool, each column is an AS event. The values in the matrix are the processed foldchange values between included exons and skipping exons read counts. 
  }
}
\details{
  This data frame is the raw output data from large-scale screening.
}

\examples{
  data(countMat)
}
\keyword{datasets}
