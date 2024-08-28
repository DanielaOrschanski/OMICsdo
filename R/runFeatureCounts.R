#' @title Run Feature Counts
#' @description Executes the feature count function from Rsubread library.
#' @param patients_dir directory that contains
#' @return FC.object featureCounts object.
#' @import Rsubread
#' @import openxlsx
#' @export
runFeatureCounts <- function(patients_dir, genomeRef) {

  list_patients <- list.dirs(path = patients_dir, full.names = FALSE, recursive = FALSE)
  cant_patients <- length(list_patients)

  list_bams <- c()
  for (p in list_patients) {
    print(p)
    file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, p))
    bam <- file_list[endsWith(file_list, "Aligned_out.bam")]

    if (length(nchar(bam)) == 0) {
      message(sprintf("There are no bam files in %s directory", p))
    } else {
      bam_file <- paste0(patients_dir, "/", p, "/trimmed/", file_list[endsWith(file_list, "Aligned_out.bam")], sep="")
      list_bams <- c(list_bams, bam_file)
    }
  }

  genome <- ifelse(genomeRef == "HG38", "hg38", "hg19")
  FC.object <- Rsubread::featureCounts(files = list_bams,
                                       annot.inbuilt = genome,
                                       juncCounts = TRUE,
                                       isPairedEnd = TRUE)
  saveRDS(FC.object, file = sprintf("%s/FeatureCount_Report.rds", patients_dir))
  message("Feature Counts' analysis has finished!")
  return(FC.object)
}
