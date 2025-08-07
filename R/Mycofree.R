####################################################################################
# PASOS PARA EL CONTROL DE CONTAMINACION EN MUESTRAS DE RNA ########################
####################################################################################

#1. Descargar el fasta del organismo contaminante desde NCBI
#Descargar el fasta de mycoplasma : https://drive.google.com/drive/u/1/folders/1-vmdp-IF8JVviz29rTgLaV94b5tjdVFS

#2. Generar index del fasta con Rsubread

#index_rsubread <- buildindex(basename = "mycoplasma_index", reference = "/media/4tb1/Daniela/Environ/Fibroblastos/Mycoplasma/Mycoplasma.fasta")

#' @title downloadBBmap
#' @description Installation of BBmap
#' @param path_for_BBMap Path of the directory where the software will be stored.
#' @return path_bbmap path of the software.
#' @export
downloadBBmap <- function(path_for_BBMap){
  path_bbmap <- sprintf("%s/bbmap", path_for_BBMap)

  if(!file.exists(path_bbmap)) {
    system(sprintf(
      'wget -O %s/BBMap_39.06.tar.gz https://sourceforge.net/projects/bbmap/files/latest/download',
      path_for_BBMap
    ))
    system(sprintf("tar -xzvf %s/BBMap_39.06.tar.gz", path_for_BBMap))
    system(sprintf("echo 'export PATH='$PATH:$%s/bbmap'' >> ~/.bashrc", path_for_BBMap))
    #system("source ~/.bashrc")
  } else {
    message("BBmap already downloaded")
  }
  return(path_bbmap)
}


#' @title runBBmap
#' @description Execution of BBmap
#' @param path_bbmap Path where the software is stored.
#' @param patient_dir Path where the patient is stored.
#' @param rate string from 0 to 1 indicating the percentage of the output fastq. For example: rate = "0.6"
#' @return list(R1_bbmap, R2_bbmap).
#' @export
#rate = "0.6"
#path_bbmap <- "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/bbmap"

runBBmap <- function(path_bbmap, patient_dir, rate) {
  id <- basename(patient_dir)
  print(id)
  file_list <- list.files(patient_dir, recursive= FALSE, full.names = TRUE)
  R1 <- file_list[endsWith(file_list, "R1.fastq.gz")]
  R2 <- file_list[endsWith(file_list, "R2.fastq.gz")]

  bbmap_dir <- sprintf("%s/AlineadoMycoplasmaRsubread", patient_dir)
  dir.create(bbmap_dir)
  R1_bbmap <- sprintf("%s/AlineadoMycoplasmaRsubread/%s_bbmap_%s_R1.fastq.gz", patient_dir, id, rate)
  R2_bbmap <- sprintf("%s/AlineadoMycoplasmaRsubread/%s_bbmap_%s_R2.fastq.gz", patient_dir, id, rate)

  if(!file.exists(R1_bbmap) | !file.exists(R2_bbmap)) {
    system(sprintf(
      "%s/reformat.sh in1=%s in2=%s out1=%s out2=%s samplerate=%s",
      path_bbmap, R1, R2, R1_bbmap, R2_bbmap, rate ))

    #system(sprintf("%s/reformat.sh in=%s out=%s samplerate=%s", path_bbmap, R2, R2_bbmap, rate))
  } else {
    message("BBmap files have already been generated")
  }

  #reformat.sh in= 5_S36_L001_R2_001.fastq.gz out= 5_S15_L001_R2_001.fastq.gz samplerate=0.6
  return(list(R1_bbmap, R2_bbmap))
}




#Alineado GRUPAL  -------------------------------------------------------------------------------

#' @title mapeoMycoplasma_Grupal
#' @description Execution of mapeoMycoplasma_Individual to every patient of the cohort.
#' @param dir_patients Path where the patients are stored.
#' @param index_rsubread Path where the index of mycoplasma generated with Rsubread is stored.
#' @param bbmap TRUE/FALSE. Set to TRUE in order to subsample the fastq files.
#' @param path_bbmap Path where the software is stored. Set to NA if BBmap is not needed.
#' @param rate string from 0 to 1 indicating the percentage of the output fastq. For example: rate = "0.6" will generate fastq files with 60% of the original size.
#' @return mapeos_mycoplasma: df with percentage of mapped reads to mycoplasma per sample.
#' @export
mapeoMycoplasma_Grupal <- function(dir_patients, index_rsubread, bbmap = FALSE, path_bbmap = NA, rate = "0.6") {

  patients <- list.dirs(dir_patients, full.names = FALSE, recursive = FALSE)
  mapeos_mycoplasma <- data.frame(ID = rep(0, length(patients)), Mapeo = rep(0,length(patients)))
  #p=2
  for (p in 1:length(patients)) {
    id <- patients[p]
    print(id)
    patient_dir <- paste(dir_patients, "/", id, sep ="")

    out <- mapeoMycoplasma_Individual(patient_dir = patient_dir,
                                      path_bbmap = path_bbmap,
                                      index_rsubread = index_rsubread,
                                      bbmap = bbmap, rate = rate)

    mapeos_mycoplasma$Mapeo[p] <- out$Mapeo
    mapeos_mycoplasma$ID[p] <- out$ID

  }
  return(mapeos_mycoplasma)
}

#--------------------------------------------------------------

#' @title mapeoMycoplasma_Individual
#' @description Execution of mapeoMycoplasma_Individual to every patient of the cohort.
#' @param patient_dir Path where the patients are stored.
#' @param index_rsubread Path where the index of mycoplasma generated with Rsubread is stored.
#' @param bbmap TRUE/FALSE. Set to TRUE in order to subsample the fastq files.
#' @param path_bbmap Path where the software is stored. Set to NA if BBmap is not needed.
#' @param rate string from 0 to 1 indicating the percentage of the output fastq. For example: rate = "0.6" will generate fastq files with 60% of the original size.
#' @return mapeo_mycoplasma: df with ID and percentage of mapped reads to mycoplasma.
#' @export
mapeoMycoplasma_Individual <- function(patient_dir, index_rsubread, bbmap = FALSE, path_bbmap = NA, rate = "0.6") {

  id <- basename(patient_dir)
  print(id)

  if(bbmap == FALSE) { #Si no se habilita el bbmap, el mycofree se hace sobre los fastq originales --> Tarda mas tiempo!
    file_list <- list.files(patient_dir, recursive= FALSE, full.names = TRUE)
    R1 <- file_list[endsWith(file_list, "R1.fastq.gz")]
    R2 <- file_list[endsWith(file_list, "R2.fastq.gz")]
    message("Mycofree va a ejecutarse sobre los files originales")
    rate <- "100"

  } else {
    out <- runBBmap(path_bbmap = path_bbmap, patient_dir = patient_dir, rate = rate)
    R1 <- out[[1]]
    R2 <- out[[2]]
    message(sprintf("Mycofree va a ejecutarse sobre los files BBmap con rate %s", rate))
  }

  #R1 <- sprintf("%s/%s/%s_R1.fastq.gz", dir_patients, patients[p], patients[p])
  #R2 <- sprintf("%s/%s/%s_R2.fastq.gz", dir_patients, patients[p], patients[p])

  path <- sprintf("%s/AlineadoMycoplasmaRsubread", patient_dir)

  library(Rsubread)
  if (file.exists(paste(path, "/", id, "_rate", rate ,"_alignRsubread.BAM.summary", sep=""))) {
    print("ya se le hizo el alineamiento")
  } else {
    print("no se le hizo el alinemiento")
    dir.create(path)
    align(
      index = index_rsubread, #Tiene que tener el nombre de todos los archivos del index, no de la carpeta nomas.
      readfile1 = R1,
      readfile2 = R2,
      input_format = "FASTQ",
      output_file = paste(path, "/",id, "_rate", rate ,"_alignRsubread.BAM", sep=""),
      type = "rna",
      nthreads = 15
    )
  }
  summary <- read.delim(paste(path, "/", id, "_rate", rate,"_alignRsubread.BAM.summary", sep=""), header = FALSE)
  total <- summary[1,2]
  mapped <- summary[3,2]
  mapeo_unico <- round(mapped/total*100, 2)

  mapeo_mycoplasma <- data.frame(ID = id, Mapeo = mapeo_unico)

  return(mapeo_mycoplasma)
}
