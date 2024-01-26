#' @title downloadFastQC
#' @description Downloads and decompresses the FastQC software
#' @return The path where the .exe file is located
#' @export
IdownloadFastQC <- function(soft_directory = sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))) {
  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
      linea_software <- grep("(?i)FastQC", softwares, ignore.case = TRUE, value = TRUE)
      FastQC <<- strsplit(linea_software, " ")[[1]][[2]]
      system2(FastQC, "--help")

      message("FastQC was already downloaded")
      #return(FastQC)
    },
    error = function(e) {
      message("FastQC download and installation is about to begin...
              Please be patient, it may take a while.")
      print(e)

      # En caso de haberlo descargado y hay algun problema con el ejecutable,
      # eliminamos y descargamos de nuevo
      if (file.exists(sprintf('%s/FastQC', soft_directory))) {
        system2("rm", sprintf('-r %s/FastQC', soft_directory))
      }

      tryCatch(
        expr = {
          # Proceso de Descarga de FastQC
          dir.create(sprintf("%s/FastQC", soft_directory))
          setwd(sprintf("%s/FastQC", soft_directory))
          URL <- ("https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip")

          system2("wget", URL, wait = TRUE, stdout = NULL, stderr = NULL)

          file <- basename(URL)
          file_dir <- file.path(sprintf("%s/FastQC", soft_directory), file)
          filedc <- substr(file, start = 0, stop = (nchar(file) - 4))

          # Descomprimo y elimino el ZIP
          system2("unzip", sprintf("%s -d %s/FastQC", file_dir, soft_directory), wait = TRUE)
          system2("rm", file_dir)

          # Definicion de la variable FastQC con el ejecutable
          FastQC <<- sprintf('%s/FastQC/FastQC/fastqc', soft_directory)
          system2(FastQC, "--help")

          # En caso de que la descarga haya sido exitosa, agregamos el path al archivo TXT
          softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))

          if(length(nchar(softwares[-grep("FastQC", softwares, ignore.case = TRUE)])) != 0 ) {
            softwares <- softwares[-grep("FastQC", softwares, ignore.case = TRUE)]
          }
          softwares_actualizado <- c(softwares, sprintf("FastQC %s", FastQC))
          # Reescribimos el archivo
          write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))


          message("FastQC download and installation completed successfully")

          return(FastQC)
        },
        error = function(e) {
          message("An error occured while performing the FastQC download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install wget unzip
------------------------------------------------------")
          print(e)
        },

        finally = {
          message("-.Message from FastQC")
        }
      )
    }

  )
}
