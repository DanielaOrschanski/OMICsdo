.onLoad <- function(libname, pkgname) {

  # HOLAAAA CAMBIE ESTO""""""""""""""""

  message("ESTOY EN EL ONLOAD")

  libPath <- dirname(system.file(package = "OMICsdo"))
  print(libPath)

  #libPath <- Sys.getenv('R_LIBS_USER')

  #Folder where the softwares will be saved
  if(!(file.exists(sprintf("%s/OMICsdoSof", libPath)))) {
    dir.create(sprintf("%s/OMICsdoSof", libPath))
  }

  if (!(file.exists(sprintf("%s/OMICsdoSof/path_to_soft.txt", libPath)))) { #Solo en instalacion
    write("",file = sprintf("%s/OMICsdoSof/path_to_soft.txt", libPath))
  }

  omicsdo_sof <<- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))

  check_packages()

  downloadFastQC()
  downloadTrimGalore()
  downloadSTAR()
  downloadArriba()
  downloadBWA()
  #downloadGATK()
  #downloadPICARD()

}

check_packages <- function() {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    remotes::install_cran("openxlsx")
  }
  library(openxlsx)

  if (!requireNamespace("readr", quietly = TRUE)) {
    remotes::install_cran("readr")
  }
  library(readr)

  if (!requireNamespace("stringr", quietly = TRUE)) {
    remotes::install_cran("stringr")
  }
  library(stringr)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    remotes::install_cran("ggplot2")
  }
  library(ggplot2)

  if (!requireNamespace("tibble", quietly = TRUE)) {
    remotes::install_cran("tibble")
  }
  library(tibble)

  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    remotes::install_cran("tidyverse")
  }
  library(tidyverse)

  if (!requireNamespace("showtext", quietly = TRUE)) {
    remotes::install_cran("showtext")
  }

  library(showtext)
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    remotes::install_cran("dplyr")
  }
  library(dplyr)

  if (!requireNamespace("ggtext", quietly = TRUE)) {
    remotes::install_cran("ggtext")
  }
  library(ggtext)

  if (!requireNamespace("readxl", quietly = TRUE)) {
    remotes::install_cran("readxl")
  }
  library(readxl)

  if (!requireNamespace("magrittr", quietly = TRUE)) {
    remotes::install_cran("magrittr")
  }
  library(magrittr)

  if (!requireNamespace("httr", quietly = TRUE)) {
    remotes::install_cran("httr")
  }
  library(httr)

  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    BiocManager::install("Rsamtools")
  }
  library(Rsamtools)

  if (!requireNamespace("viridis", quietly = TRUE)) {
    install.packages("viridis")
  }
  library(viridis)

  if (!requireNamespace("reshape2", quietly = TRUE)) {
    install.packages("reshape2")
  }
  library(reshape2)

  #Para usar MIXTURE
  if (!requireNamespace("nnls", quietly = TRUE)) {
    install.packages("nnls")
  }
  library(nnls)

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    BiocManager::install("ComplexHeatmap", force = TRUE)
  }
  library(ComplexHeatmap)

  if (!requireNamespace("MIXTURE", quietly = TRUE)) {
    install_github("elmerfer/MIXTURE")
  }
  library(MIXTURE)

  if (!requireNamespace("Rsubread", quietly = TRUE)) {
    BiocManager::install("Rsubread")
  }
  library(Rsubread)

  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery")
  }
  library(GEOquery)

  message("The R packages required have been successfully installed")
}

#' @title downloadFastQC
#' @description Downloads and decompresses the FastQC software
#' @return The path where the .exe file is located
downloadFastQC <- function() {
  soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo")) ))
      linea_software <- grep("(?i)FastQC", softwares, ignore.case = TRUE, value = TRUE)

      if(length(nchar(linea_software))==0) {
        stop()
      } else {
        FastQC <<- strsplit(linea_software, " ")[[1]][[2]]
        system2(FastQC, "--help")
        message("FastQC was already downloaded")
        return(FastQC)
      }

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
          softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

          if(length(nchar(softwares[-grep("FastQC", softwares, ignore.case = TRUE)])) != 0 ) {
            softwares <- softwares[-grep("FastQC", softwares, ignore.case = TRUE)]
          }
          softwares_actualizado <- c(softwares, sprintf("FastQC %s", FastQC))
          # Reescribimos el archivo
          write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))


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

######################################################################
######################################################################
#' @title downloadTrimGalore
#' @description Downloads and decompresses the TrimGalore software
#' @return The path where the .exe file is located
downloadTrimGalore <- function() {

  soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      linea_software <- grep("(?i)TrimGalore", softwares, ignore.case = TRUE, value = TRUE)

      if(length(nchar(linea_software))==0) {
        stop()
      } else {
        TrimGalore <<- strsplit(linea_software, " ")[[1]][[2]]
        system2(TrimGalore, "--help")
        message("TrimGalore was already downloaded")
        return(TrimGalore)
      }

    },
    error = function(e) {
      # En caso de haberlo descargado y hay algun problema con el ejecutable, eliminamos y descargamos de nuevo
      if (file.exists(sprintf('%s/TrimGalore', soft_directory))) {
        system2("rm", sprintf('-r %s/TrimGalore', soft_directory))
        print("There is a problem with the TrimGalore exe file. It will be removed and download again")
      }

      message("TrimGalore download and installation is about to begin...
              Please be patient, it may take a while.")
      print(e)

      tryCatch(
        expr = {
          # Primero hay que verificar que este descargado FastQC y Cutadapt.
          downloadFastQC()

          # Proceso de Descarga de TrimGalore
          dir.create(sprintf("%s/TrimGalore", soft_directory))
          URL <- "https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz"
          setwd(sprintf("%s/TrimGalore", soft_directory))
          system2("wget" ,sprintf("https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -O trim_galore.tar.gz"), wait = TRUE, stdout = NULL, stderr = NULL)
          system2("tar", "xvzf trim_galore.tar.gz")

          TrimGalore <<- sprintf("%s/TrimGalore/TrimGalore-0.6.10/trim_galore", soft_directory)
          system2(TrimGalore, "--help")

          # En caso de que la descarga haya sido exitosa, agregamos el path al archivo TXT
          softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
          if(length(nchar(softwares[-grep("TrimGalore", softwares, ignore.case = TRUE)])) != 0 ) {
            softwares <- softwares[-grep("TrimGalore", softwares, ignore.case = TRUE)]
          }
          softwares_actualizado <- c(softwares, sprintf("TrimGalore %s", TrimGalore))
          write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

          return(TrimGalore)
        },

        error = function(e) {
          message("An error occured while performing the TrimGalore download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install cutadapt wget tar
------------------------------------------------------")
          print(e)
        },

        finally = {
          message("-.Message from TrimGalore")
        }
      )
    },
    finally = {
      message("TrimGalore download and installation completed successfully")
    }
  )
}
##############################################################################
##############################################################################
#' @title downloadSTAR
#' @description Downloads and decompresses the GATK software
#' @return The path where the .exe file is located
downloadSTAR <- function() {
  soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))

  tryCatch(
    expr = {

      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      #Para que sea insensible a las mayusculas o minusculas se pone el (?i)
      linea_software <- grep("(?i)STAR", softwares, ignore.case = TRUE, value = TRUE)

      if(length(nchar(linea_software))== 0) {
        stop()
      } else {
        STAR <<- strsplit(linea_software, " ")[[1]][[2]]
        system2(STAR, "--help")
        return(STAR)
      }

    },
    error = function(e) {

      message("ESTOY EN EL ERROR DEL STAR")
      # En caso de haberlo descargado y hay algun problema con el ejecutable,
      # eliminamos y descargamos de nuevo
      if (file.exists(sprintf('%s/STAR', soft_directory))) {
        system2("rm", sprintf('-r %s/STAR', soft_directory))
        print("There is a problem with the STAR exe file. It will be removed and download again")
      }

      message("STAR download and installation is about to begin...
              Please be patient, it may take a while.")
      print(e)

      tryCatch(
        expr = {
          # Proceso de Descarga de STAR
          dir.create(sprintf("%s/STAR", soft_directory))
          URL <- "https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz"
          version <- sub("\\.tar\\.gz$", "", basename(URL))
          dir <- sprintf("%s/STAR", soft_directory)
          setwd(dir)
          # Get latest STAR source from releases
          system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
          system2("tar", sprintf("-xzf %s.tar.gz", version))
          # Compile
          setwd(sprintf("%s/STAR-%s/source", dir, version))
          system2("make", "STAR")

          # Definimos la variable de STAR con su ejecutable
          STAR <<- sprintf("%s/STAR/STAR-2.7.11a/bin/Linux_x86_64_static/STAR", soft_directory)
          system2(STAR, "--help")

          # En caso de que la descarga haya sido exitosa, agregamos el path al archivo TXT
          softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
          softwares_actualizado <- c(softwares, sprintf("STAR %s", STAR))
          write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

          return(STAR)
        },

        error = function(e) {
          message("An error occured while performing the STAR download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
    sudo apt install wget tar
------------------------------------------------------")
          print(e)
        },

        finally = {
          message("-.Message from STAR")
        }
      )
    },
    finally = {
      message("STAR download and installation completed successfully")
    }
  )
}


####################################################################
######################################################################

#' @title downloadArriba
#' @description Downloads and decompresses the GATK software
#' @return The path where the .exe file is located

downloadArriba <- function() {
  soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))

  tryCatch(
    expr = {
      softwares <- readLines(sprintf("%s/path_to_soft.txt", soft_directory))
      linea_software <- grep("(?i)Arriba", softwares, ignore.case = TRUE, value = TRUE)

      if(length(nchar(linea_software)) == 0) {
        stop()
      } else {
        ARRIBA <<- strsplit(linea_software, " ")[[1]][[2]]
        system2(ARRIBA, "-help")
        return(ARRIBA)
      }

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
          softwares <- readLines(sprintf("%s/path_to_soft.txt",soft_directory))
          if (TRUE %in% grepl("ARRIBA", softwares, ignore.case = TRUE)) {
            # En caso de haber dado error y se descargo de nuevo, tenemos que eliminar la linea del
            # software anterior.
            softwares <- softwares[-grep("Arriba", softwares, ignore.case = TRUE)]
          }
          # Agregamos la nueva linea con el software
          softwares_actualizado <- c(softwares, sprintf("ARRIBA %s", ARRIBA))
          write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

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




################################################################
# CNVS ##################################3
########################################################


# BWA
#' @title downloadBWA
#' @description Downloads and decompresses the BWA software
#' @return The path where the .exe file is
downloadBWA <- function() {
  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  print(omicsdo_sof)
  #omicsdo_sof <- sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))

  tryCatch(
    expr = {
      system(sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof))
    },
    error = function(e) {
      message("Installation of BWA will now begin. Check if all the required packages are downloaded.")
      print(e)

      dir.create(sprintf("%s/BWA", omicsdo_sof))

      bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
      bwa_dir1 <- file.path(omicsdo_sof, "BWA")
      system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", omicsdo_sof))
      bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
      bwa_dir2 <- file.path(omicsdo_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
      system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)

      bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
      system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      BWA <- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("BWA %s", BWA))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },
    warning = function(w) {
      message("Installation of BWA will now begin. Check if all the required packages are downloaded.")

      dir.create(sprintf("%s/BWA", omicsdo_sof))

      bwa_url1 <- "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm"
      bwa_dir1 <- file.path(omicsdo_sof, "BWA")
      system2("wget", args = c(bwa_url1, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      dir.create(sprintf("%s/BWA/bwa-0.7.17-lp154.6.1.src", omicsdo_sof))
      bwa_url2 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm"
      bwa_dir2 <- file.path(omicsdo_sof, "BWA/bwa-0.7.17-lp154.6.1.src")
      system2("wget", args = c(bwa_url2, "-P", bwa_dir2), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.src.rpm | cpio -D %s -idmv", bwa_dir2, bwa_dir2), wait = TRUE)

      bwa_url3 <-"https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm"
      system2("wget", args = c(bwa_url3, "-P", bwa_dir1), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("rpm2cpio", sprintf("%s/bwa-0.7.17-lp154.6.1.i586.rpm | cpio -D %s -idmv", bwa_dir1, bwa_dir1), wait = TRUE)

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      BWA <- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("BWA %s", BWA))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },

    finally = {
      BWA <<- sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof)
      message("-.Message from BWA")
    }
  )

  return(sprintf('%s/BWA/usr/bin/bwa', omicsdo_sof))
}

# GATK
#' @title downloadGATK
#' @description Downloads and decompresses the GATK software
#' @return The path where the .exe file is located
downloadGATK <- function() {
  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  #omicsdo_sof <- sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))
  tryCatch(
    expr = {
      system(sprintf('java -jar %s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof))
    },
    error = function(e) {
      message("Installation of GATK will now begin. Check if all the required packages are downloaded.")
      print(e)
      dir.create(sprintf("%s/GATK", omicsdo_sof))
      URL <- 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip'
      system2("wget", args = c(URL, "-P", paste(omicsdo_sof, "/GATK", sep="")), wait = TRUE, stdout = NULL, stderr = NULL)
      file <- basename(URL)
      file_dir <- paste(omicsdo_sof, "/GATK/", file , sep = "")
      filedc <- substr(file, start = 0, stop = (nchar(file) - 4))
      unzip(zipfile = file_dir, exdir = paste(omicsdo_sof, "/GATK", sep=""))

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      GATK <- sprintf('%s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("GATK %s", GATK))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },
    warning = function(w) {
      message("Installation of GATK will now begin. Check if all the required packages are downloaded.")
      dir.create(sprintf("%s/GATK", omicsdo_sof))
      URL <- 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip'
      system2("wget", args = c(URL, "-P", paste(omicsdo_sof, "/GATK", sep="")), wait = TRUE, stdout = NULL, stderr = NULL)
      file <- basename(URL)
      file_dir <- paste(omicsdo_sof, "/GATK/", file , sep = "")
      filedc <- substr(file, start = 0, stop = (nchar(file) - 4))
      unzip(zipfile = file_dir, exdir = paste(omicsdo_sof, "/GATK", sep=""))

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      GATK <- sprintf('%s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("GATK %s", GATK))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },
    finally = {
      GATK <<- sprintf('%s/GATK/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar', omicsdo_sof)
      message("-.Message from GATK")
    }
  )

  return(sprintf('%s/GATK/%s/gatk-package-4.3.0.0-local.jar', omicsdo_sof, filedc))
}


# PICARD
#' @title downloadPICARD
#' @description Downloads and decompresses the PICARD software
#' @return The path where the .exe file is located
downloadPICARD <- function() {

  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))

  tryCatch(
    expr = {
      system(sprintf('java -jar %s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof))
    },
    error = function(e) {
      message("Installation of PICARD will now begin. Check if all the required packages are downloaded.")
      print(e)
      dir.create(sprintf("%s/PICARD", omicsdo_sof))
      URL <- "https://github.com/broadinstitute/picard/archive/refs/tags/2.27.5.tar.gz"
      system2("wget", args = c(URL, "-P", paste(omicsdo_sof, "/PICARD", sep="")), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("gzip" , sprintf("-d %s/PICARD/2.27.5.tar.gz", omicsdo_sof))
      system2("tar" , sprintf("-xvf %s/PICARD/2.27.5.tar -C %s/PICARD", omicsdo_sof, omicsdo_sof))
      file.remove(sprintf("%s/PICARD/2.27.5.tar", omicsdo_sof))

      URL2 <- "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar"
      dir2 <- sprintf("%s/PICARD/picard-2.27.5", omicsdo_sof)
      system2("wget", args = c(URL2, "-P", dir2), wait = TRUE, stdout = NULL, stderr = NULL)

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      PICARD <- sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("PICARD %s", PICARD))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },
    warning = function(w) {
      message("Installation of PICARD will now begin. Check if all the required packages are downloaded.")
      dir.create(sprintf("%s/PICARD", omicsdo_sof))
      URL <- "https://github.com/broadinstitute/picard/archive/refs/tags/2.27.5.tar.gz"
      system2("wget", args = c(URL, "-P", paste(omicsdo_sof, "/PICARD", sep="")), wait = TRUE, stdout = NULL, stderr = NULL)
      system2("gzip" , sprintf("-d %s/PICARD/2.27.5.tar.gz", omicsdo_sof))
      system2("tar" , sprintf("-xvf %s/PICARD/2.27.5.tar -C %s/PICARD", omicsdo_sof, omicsdo_sof))
      file.remove(sprintf("%s/PICARD/2.27.5.tar", omicsdo_sof))

      URL2 <- "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar"
      dir2 <- sprintf("%s/PICARD/picard-2.27.5", omicsdo_sof)
      system2("wget", args = c(URL2, "-P", dir2), wait = TRUE, stdout = NULL, stderr = NULL)

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      PICARD <- sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("PICARD %s", PICARD))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },

    finally = {
      PICARD <<- sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof)
      message("-.Message from PICARD")
    }
  )

  return(sprintf('%s/PICARD/picard-2.27.5/picard.jar', omicsdo_sof))
}


# SAMTOOLS
#' @title downloadSamtools
#' @description Downloads and decompresses the Samtools software
#' @return The path where the .exe file is located
downloadSamtools <- function() {

  tryCatch(
    {
      system2(sprintf("%s/Samtools/samtools-1.16.1/samtools", omicsdo_sof))
    },
    error = function(e) {
      message("The installation of Samtools will begin now. Check if all required packages are downloaded.")
      print(e)

      dir.create(sprintf("%s/Samtools", omicsdo_sof))
      dir <- sprintf("%s/Samtools", omicsdo_sof)
      URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
      system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)

      system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", omicsdo_sof))
      system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", omicsdo_sof, list.files(sprintf("%s/Samtools", omicsdo_sof)), dir)))

      file.remove(sprintf("%s/samtools-1.16.1.tar", dir))

      samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
      system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
      system(paste("cd", shQuote(samtools_dir), "&& make"))

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      Samtools <- sprintf("%s/Samtools/samtools-1.16.1/samtools", omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("Samtools %s", Samtools))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },
    warning = function(w) {
      message("The installation of Samtools will begin now. Check if all required packages are downloaded.")

      dir.create(sprintf("%s/Samtools", omicsdo_sof))
      dir <- sprintf("%s/Samtools", omicsdo_sof)
      URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
      system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)

      system2("bzip2", sprintf("-d %s/Samtools/samtools-1.16.1.tar.bz2", omicsdo_sof))
      system2("tar", c("-xvf", sprintf("%s/Samtools/%s -C %s", omicsdo_sof, list.files(sprintf("%s/Samtools", omicsdo_sof)), dir)))

      file.remove(sprintf("%s/samtools-1.16.1.tar", dir))

      samtools_dir <- sprintf("%s/samtools-1.16.1", dir)
      system(paste("cd", shQuote(samtools_dir), "&& ./configure"))
      system(paste("cd", shQuote(samtools_dir), "&& make"))

      #Escribo el path en el txt
      softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
      Samtools <- sprintf("%s/Samtools/samtools-1.16.1/samtools", omicsdo_sof)
      softwares_actualizado <- c(softwares, sprintf("Samtools %s", Samtools))
      write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    },

    finally = {
      message("-.Message from Samtools")
      Samtools <<- sprintf("%s/Samtools/samtools-1.16.1/samtools", omicsdo_sof)
    }
  )

  return(sprintf("%s/Samtools/samtools-1.16.1/samtools", omicsdo_sof))
}



