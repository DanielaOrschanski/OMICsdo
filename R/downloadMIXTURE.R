#' @title downloadHG38
#' @description Downloads the FASTA and the annotation of the genome reference version HG38.
#' @return Paths of FASTA and GTF(annotation) from the genome reference.
#' @import GEOquery
#' @export

downloadMIXTURE <- function() {

  install.packages("nnls")
  library(nnls)
  BiocManager::install("ComplexHeatmap", force = TRUE)
  install_github("elmerfer/MIXTURE")
  library(MIXTURE)
  data(LM22)
  data("TIL10")
  LM22 <- LM22
  TIL10 <- TIL10
}



