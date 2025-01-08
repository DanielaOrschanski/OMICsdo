getwd()
setwd("/media/16TBDisk")

setwd("/media/16TBDisk/Daniela")

# Instala el paquete googledrive si no lo has instalado aún
if (!requireNamespace("googledrive", quietly = TRUE)) {
  install.packages("googledrive")
}
library(googledrive)

drive_auth()

# Para descargar un solo archivo ---------------------------------------------------------------------

library(googledrive)

url <- "https://drive.google.com/file/d/1-vxo2G5p62KtlPTJwg34Rf4IYZlQbGCa/view?usp=drive_link"

file_id <- basename(dirname(url))

local_path <- "/media/4tb1/Daniela/CNVs/MuestrasCNVs/Axen/FNM-020-3gb.rawdata"
download_url <- paste0("https://drive.google.com/uc?id=", file_id)

drive_download(as_id(file_id), path = local_path, overwrite = TRUE)




url <- "https://drive.google.com/file/d/1B3xIadGGo5WkqG-IdDHYyW7N22zbejPG/view?usp=drive_link"
file_id <- basename(dirname(url))
local_path <- "/media/16TBDisk/Daniela/CNVs/MuestrasCNVs/37513/OMICsdo/37513_recal.bai"
download_url <- paste0("https://drive.google.com/uc?id=", file_id)
drive_download(as_id(file_id), path = local_path, overwrite = TRUE)


url <- "https://drive.google.com/file/d/1ojJ3L6iXcSU-NimFFHaEekbWkcMUkmTV/view?usp=drive_link"
file_id <- basename(dirname(url))
local_path <- "/media/16TBDisk/Daniela/CNVs/MuestrasCNVs/37821/37821_recal.bai"
download_url <- paste0("https://drive.google.com/uc?id=", file_id)
drive_download(as_id(file_id), path = local_path, overwrite = TRUE)


#Para descargar varios archivos de una carpeta en Drive ------------------------

# URL de la carpeta
folder_url <- "https://drive.google.com/drive/folders/1z7bh5cB9Iu3oYYoCNvqtJIJixxXNPHnS?usp=drive_link"
#Tiene que tener la barra / al final del path!!!!!!!!
local_folder <- "/media/4tb2/Daniela/Biota/Nuevas11_04/"

downloadFolderDrive(folder_url,
                    local_folder)
library(googledrive)

downloadFolderDrive <- function(folder_url,
                                local_folder
                                ) {

  folder_id <- sub(".*/folders/(.*)\\?usp=drive_link", "\\1", folder_url)
  folder_files <- as.data.frame(drive_ls(as_id(folder_id)))

  # Iterar sobre la lista de archivos y descargarlos

  for (i in 1:nrow(folder_files)) {
    file_id <- folder_files[i,"id"]
    file_name <- folder_files[i,"name"]

    #Crea una carpeta donde se van a guardar las muestras
    id <- strsplit(file_name, split="_")[[1]][1]
    dir.create(paste(local_folder, id, sep=""))

    gz <- ifelse(grepl("gz", file_name), ".gz", "")
    R <- ifelse(grepl("R1", file_name), "R1", "R2")
    local_path <- paste(local_folder, id, "/", id, "_S04_L001_", R, "_001.fastq", gz, sep="")
    #setwd(local_path)
    drive_download(file_id, path = local_path, overwrite = TRUE)
  }

}


#Para descargar muestras de CNVs de fleni -----------------------------------------------------
library(googledrive)
folder_url <- "https://drive.google.com/drive/folders/1eKG6m7cGAVOuLqKZjBn8G3jpau76aHBt?usp=drive_link"
#Tiene que tener la barra / al final del path!!!!!!!!
local_folder <- "/media/4tb2/Daniela/Mama TNF HER2/Muestras/"

downloadFolderDrive(folder_url, local_folder)

downloadFolderDrive <- function(folder_url,
                                local_folder
) {

  folder_id <- sub(".*/folders/(.*)\\?usp=drive_link", "\\1", folder_url)
  folder_files <- as.data.frame(drive_ls(as_id(folder_id)))

  # Iterar sobre la lista de archivos y descargarlos

  for (i in 1:nrow(folder_files)) {
    file_id <- folder_files[i,"id"]
    file_name <- folder_files[i,"name"]

    #Crea una carpeta donde se van a guardar las muestras
    patron <- ".*-([^-_]+).*"
    id <- sub(patron, "\\1", file_name)
    dir.create(paste(local_folder, id, sep=""))
    dir.create(paste(local_folder, id, "/OMICsdo",  sep=""))

    gz <- ifelse(grepl("gz", file_name), ".gz", "")

    #Se fija si es R1, R2, bam o bai:
    if (grepl("R1", file_name)) {
      local_path <- paste(local_folder, id, "/", id, "_R1.fastq", gz, sep="")
    } else if (grepl("R2", file_name)) {
      local_path <- paste(local_folder, id, "/", id, "_R2.fastq", gz, sep="")
    } else if (grepl("bam", file_name)) {
      if(grepl("bai", file_name)) {
        local_path <- paste(local_folder, id, "/OMICsdo", id, "_sortedR1R2.bam.bai", sep="")
      }
      local_path <- paste(local_folder, id, "/OMICsdo", id, "_sortedR1R2.bam", sep="")
    }

    drive_download(file_id, path = local_path, overwrite = TRUE)
  }

}

# Subir archivos -------------------------
list_pat <- c("BM22-46332", "BM22-46248", "BM22-46334", "BM22-46432" )
for(i in list_pat)  {
  drive_upload(sprintf("/home/biomolecular/DATA/NGS/Pacientes/%s_MODApy/%sMODApy_realigned_reads_recal.bam", i, i), name = sprintf("%s_sortedR1R2.bam", i))
  drive_upload(sprintf("/home/biomolecular/DATA/NGS/Pacientes/%s_MODApy/%sMODApy_realigned_reads_recal.bai", i, i), name = sprintf("%s_sortedR1R2.bam.bai", i))
}

i="32"
drive_upload(sprintf("/media/16TBDisk/Daniela/Biota/Muestras2daTanda/5/concat5_%s/concat5_%s_R1.fastq", i, i), name = sprintf("concat5_%s_R1.fastq", i))
drive_upload(sprintf("/media/16TBDisk/Daniela/Biota/Muestras2daTanda/5/concat5_%s/concat5_%s_R2.fastq", i, i), name = sprintf("concat5_%s_R2.fastq", i))


#######################################3
library(googledrive)
folder_url <- "https://drive.google.com/drive/folders/1eKG6m7cGAVOuLqKZjBn8G3jpau76aHBt?usp=drive_link"
#Tiene que tener la barra / al final del path!!!!!!!!
local_folder <- "/media/4tb2/Daniela/Mama TNF HER2/Muestras/"

downloadFolderDrive(folder_url, local_folder)

downloadFolderDrive <- function(folder_url,
                                local_folder
) {

  folder_id <- sub(".*/folders/(.*)\\?usp=drive_link", "\\1", folder_url)
  folder_files <- as.data.frame(drive_ls(as_id(folder_id)))

  # Iterar sobre la lista de archivos y descargarlos

  for (i in 1:nrow(folder_files)) {
    file_id <- folder_files[i,"id"]
    file_name <- folder_files[i,"name"]

    #Crea una carpeta donde se van a guardar las muestras
    #patron <- ".*-([^-_]+).*"
    #id <- sub(patron, "_S", file_name)
    id <- strsplit(file_name, split= "_S")[[1]][1]
    dir.create(paste(local_folder, id, sep=""))

    gz <- ifelse(grepl("gz", file_name), ".gz", "")

    local_path <- sprintf("%s%s/%s", local_folder, id, file_name)
    drive_download(file_id, path = local_path, overwrite = TRUE)
  }

}



