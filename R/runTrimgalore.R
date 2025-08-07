#' @title RunTrimgalore
#' @description Corre la funcion TrimGalore para eliminar los adapters y las lecturas (y bases) de baja calidad.
#' @param patient_dir Path donde se encuentra el archivo R1 de formato fasta o fastq.
#' @param trim_quality is the minimum value of the quality of each base within the sequence that will pass the filter.
#' @return path of the trimmed folder which contained the trimmed files and was created inside the patient folder.
#' @export
runTrimgalore <- function(patient_dir,  trim_quality = 30, soft_directory) {

  # Chequeamos que esta descargado TrimGalore. En caso de no estarlo, lo descarga
  TrimGalore <- downloadTrimGalore(soft_directory)
  patient_id <- basename(patient_dir)
  file_list <- list.files(patient_dir)

  gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
  fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))], sep="")
  fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))], sep="")
  # Se fija si los archivos de entrada son de formato .gz
  gziped <- ifelse(stringr::str_detect(fileR1,".gz"),"--gzip","--dont_gzip")

  if ((length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0) | (length(nchar(file_list[endsWith(file_list, "R2.fastq.gz")])) == 0)) {
    message("There are no fastq files in this directory")
  }

  #Se fija si está hecho el trimmeado antes:
  trim_dir <- list.files(paste0(patient_dir, "/trimmed"))
  if (!(length(nchar(file_list[endsWith(trim_dir, "val_1.fq.gz")])) == 0)) {
    #ifelse(gziped == "--gzip",
    #       file_trim <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1.fq.gz")], sep=""),
    #       file_trim <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1.fq")], sep=""))

    message("You have already trimmed this sample")

    return(sprintf("%s/trimmed", patient_dir))
  } else {
    dir.create(sprintf("%s/trimmed", patient_dir))
    outdir <- sprintf("%s/trimmed", patient_dir)

    # Trimeado por consola. Guarda el tiempo que tardo en ejecutarse
    # PAIRED - END
    t1 <- system.time(system2(command = TrimGalore,
                              args =c("--paired",
                                      gziped,
                                      sprintf("-q %s", trim_quality),
                                      paste0("--output_dir ", outdir),
                                      fileR1,
                                      fileR2,
                                      paste0("--basename ", basename(patient_dir)))))

    #stderr = file.path(outdir, "ErrLog.txt")))

    #stderr = file.path(outdir, "ErrLog.txt")))

    # Calculo el tamaño de los archivos de entrada sin trimmear
    of.size1 <- file.info(fileR1)$size
    of.size2 <- file.info(fileR2)$size

    # Si el formato era .gz, el output se llamara tambien como .gz, porque TrimGalore asi lo genera
    if (gziped == "--gzip") {
      ofile <- paste0(outdir, "/", basename(patient_dir),"_val_1.fq.gz")
    } else {
      ofile <- paste0(outdir, "/", basename(patient_dir),"_val_1.fq")
    }

    # Tamanos del archivo de salida
    ot.size1 <- file.info(ofile)$size
    ot.size2 <- file.info(stringr::str_replace_all(ofile,"val_1.","val_2."))$size

    # No entiendo porque hace de return esta asignacion de variables
    #return(c(ifile = fileR1, tfile = ofile , Time = t1, isize1 = of.size1, isize2 = of.size2,
    #         tsize1 = ot.size1, tsize2 = ot.size2))
    message("TrimGalore's analysis has finished!")
    return(dirname(ofile))
  }

}
