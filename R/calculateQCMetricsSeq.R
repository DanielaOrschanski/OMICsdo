#' @title Calculate QC metrics from sequencing - FastQC
#' @description Executes tkjvnfkd.
#' @param patients_dir Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder.
#' @param trimmed kjvdfk
#' @return scores_qc df with metrics from FastQC
#' @export

calculateQCMetricsSeq <- function(patients_dir, trimmed = FALSE){
  library(qckitfastq)
  scores_qc <- data.frame("Sample"= c(), "Q_Mean_R1" = c(), "%_>=Q30_R1"= c(),
                          "Q_Mean_R2" = c(), "%_>=Q30_R2"= c())
  i=1
  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)

  for (p in dir_list) {
    #p <- dir_list[[i]]
    print(p)
    if  (trimmed == TRUE) {
      file_list <- list.files(sprintf("%s/trimmed/", p))
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "val_1.fq.gz")])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("val_1.fq%s", gzip))], sep="")
      fileR2 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("val_2.fq%s", gzip))], sep="")
    } else {
      file_list <- list.files(p)
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/", file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))], sep="")
      fileR2 <- paste0(p, "/", file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))], sep="")
    }

    if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
      stop("There are no fastq files in this directory")
    }

    QC <- qual_score_per_read(fileR1)
    mean <- mean(QC$mu_per_read)
    Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    scores_qc[i, "Sample"] <- basename(p)
    scores_qc[i, "Q_Mean_R1"] <- mean
    scores_qc[i, "%_>=Q30_R1"] <- Q30

    QC <- qual_score_per_read(fileR2)
    mean <- mean(QC$mu_per_read)
    Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    scores_qc[i, "Q_Mean_R2"] <- mean
    scores_qc[i, "%_>=Q30_R2"] <- Q30

    i= i+1
  }
  #saveRDS(scores_qc, file= sprintf("%s/scores_QC.rds", patients_dir))
  trimeado <- ifelse(trimmed == TRUE, "trimmed", "")
  write.xlsx(scores_qc, file= paste(patients_dir, "/MetricasQC_", trimeado,".xlsx", sep =""), rowNames= TRUE)

  return(scores_qc)

}


