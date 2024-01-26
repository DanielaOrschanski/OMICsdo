#' @title downloadArriba
#' @description Downloads and decompresses the GATK software
#' @return The path where the .exe file is located
#' @export
IdownloadARRIBA <- function() {
  soft_directory <- sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))

  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
      linea_software <- grep("(?i)Arriba", softwares, ignore.case = TRUE, value = TRUE)
      ARRIBA <<- strsplit(linea_software, " ")[[1]][[2]]

      system2(ARRIBA, "-help")
      return(ARRIBA)
    },
    error = function(e) {
      # En caso de haberlo descargado y hay algun problema con el ejecutable,
      # eliminamos y descargamos de nuevo
      if (file.exists(sprintf('%s/Arriba', soft_directory))) {
        system2("rm", sprintf('-r %s/Arriba', soft_directory))
        print("There is a problem with the Arriba .exe file. It will be removed and download again")
      }

      message("Arriba download and installation is about to begin...
              Please be patient, it may take a while.")
      print(e)

      tryCatch(
        expr = {
          # Proceso de Descarga de Arriba
          dir.create(sprintf("%s/Arriba", soft_directory))
          URL <- "https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz"
          setwd(sprintf("%s/Arriba", soft_directory))

          system2("wget" , URL, wait = TRUE, stdout = NULL, stderr = NULL)
          system2("tar", "-xzf arriba_v2.4.0.tar.gz", wait = TRUE, stdout = NULL, stderr = NULL)
          setwd(sprintf("%s/Arriba/arriba_v2.4.0", soft_directory))
          system2("make", wait = TRUE, stdout = NULL, stderr = NULL)


          ARRIBA <<- sprintf("%s/Arriba/arriba_v2.4.0/arriba", soft_directory)
          system2(ARRIBA, "-help")

          # En caso de ya existir el archivo, solamente agregamos el proximo path
          softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
          if (TRUE %in% grepl("ARRIBA", softwares, ignore.case = TRUE)) {
            # En caso de haber dado error y se descargo de nuevo, tenemos que eliminar la linea del
            # software anterior.
            softwares <- softwares[-grep("Arriba", softwares, ignore.case = TRUE)]
          }
          # Agregamos la nueva linea con el software
          softwares_actualizado <- c(softwares, sprintf("ARRIBA %s", ARRIBA))
          write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))

          return(ARRIBA)
        },

        error = function(e) {
          message("An error occured while performing the Arriba download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install wget tar
------------------------------------------------------")
          print(e)
        },

        finally = {
          message("-.Message from Arriba")
        }
      )
    },
    finally = {
      message("Arriba download and installation completed successfully")
    }
  )
}

