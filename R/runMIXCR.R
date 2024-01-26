#' @title runMIXCR
#' @description Runs the MIXCR software
#' @param patient_dir the path of the folder´s patient yoy want to analyze. It can be the one that contains the R1 and R2 original or the trimmed files.
#' @return The path where the executable file is located
#' @export

runMIXCR <- function(patient_dir) {

  file_list <- list.files(patient_dir)

  #Para que se pueda poner como entrada la carpeta de los trimmeados o la carpeta original:
  if  (startsWith(basename(patient_dir), "trimmed")) {
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "val_1.fq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("val_1.fq%s", gzip))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("val_2.fq%s", gzip))], sep="")
    patient_id <- basename(dirname(patient_dir))
  } else {
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))], sep="")
    patient_id <- basename(patient_dir)
  }

  if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
    stop("There are no fastq files in this directory")
  }

  #Evita repetir el analisis si ya fue hecho
  if (!(length(nchar(file_list[endsWith(file_list, sprintf("MIXCR_%s", patient_id))])) == 0)) {
    message("The MIXCR's analysis for this sample has already been done.")
    return(paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("MIXCR_%s", patient_id))], sep=""))
  }

  dir.create(sprintf("%s/MIXCR", patient_dir))

  #Probar esta:
  #t1 <- (system2(command =  MIXCR,
  #                          args = c("analyze rna-seq",
  #                                   "--species hsa",
  #                                   fileR1,
  #                                   fileR2,
  #                                   sprintf("%s/MIXCR/MIXCR_%s", patient_dir, patient_id)
  #                                   )




  system2(command = MIXCR,
          args = c("analyze",
            "takara-human-bcr-full-length",
                   fileR1,
                   fileR2,
                   sprintf("%s/MIXCR/%s", patient_dir, patient_id)
                   ))

  message("MIXCR's analysis has finished!")
  return(sprintf("%s/MIXCR", patient_dir))

}
