#' @title Load bed file 
#' @description Load a dataset with information about the 37 genes within the mitochondrial genome.
#' @examples  data("bedfileMito")
#' @export

bedfileHG38 <- function() {
  #load(paste(Sys.getenv('R_LIBS_USER'), "/OMICSdo/data/bedfileHG38.RData", sep=""))
  load("/media/16TBDisk/Daniela/CNVs/bedHG38.RData")
}


