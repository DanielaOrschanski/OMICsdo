#' @title Obtain statistics from fusions.
#' @description Generates 2 dataframes: 1. Rbind from every fusion report of the samples and 2. Stats from all samples indicating group.
#' @param patients_dir path to the folder that contains one folder for each sample.
#' @param Metadata dataframe that contains information about all the samples. One sample per row.
#' @param group must be a name of one of the columns of the metadata
#' @return Todos_FusionReport and Stats_Fusions: 2 dataframes.
#' @export
#' @import openxlsx
#' @import readxl
#' @examples fusionStats(patients_dir, Metadata, group = "group")

fusionStats <- function(patients_dir, Metadata, group) {

  ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  length(ids)
  #ids <- ids[-1]

  #MetadataSRA <- read.table("/media/4tb1/Daniela/Environ/Fusiones/SraRunTable.txt", header = TRUE, sep = ",")
  #colnames(MetadataSRA)[which(colnames(MetadataSRA) == "metastasis")] <- "MTT"


  Todos_FusionReport <- data.frame()
  Stats_Fusions <- data.frame()
  k=1
  for (id in ids) {
    i <- basename(id)
    print(i)
    fusions_file <- sprintf("%s/trimmed/%s_FusionReport.xlsx", id, i)
    FusionReport <- read_excel(fusions_file)

    Grupo <- Metadata[which(Metadata$ID == i), group]

    #FusionReport$MTT <- Metadata[which(Metadata$ID == i), "MTT"]
    FusionReport$Grupo <- Grupo
    FusionReport$ID <-i

    FusionReport_ID <- cbind(FusionReport$ID, FusionReport[,1:(ncol(FusionReport)-1)])
    colnames(FusionReport)[1] <- "ID"

    if(nrow(Todos_FusionReport) == 0) {
      Todos_FusionReport <- FusionReport_ID
    } else {
      Todos_FusionReport <- rbind(Todos_FusionReport, FusionReport_ID)
    }


    #Estadisticas de los reportes -----------------------------------

    Stats_Fusions[k, "ID"] <- i
    Stats_Fusions[k, "Cantidad_Fusiones"] <- nrow(FusionReport)

    if(!(nrow(FusionReport) == 0)) {
      contador_high <- 0
      contador_medium <- 0
      contador_low <- 0
      for (j in 1:nrow(FusionReport)) {
        if(FusionReport$confidence[j] == "high"){
          contador_high <- contador_high + 1
        } else if (FusionReport$confidence[j] == "medium"){
          contador_medium <- contador_medium + 1
        } else {
          contador_low <- contador_low + 1
        }
      }

      Stats_Fusions[k, "Fusiones_conf_H"] <- contador_high
      Stats_Fusions[k, "Fusiones_conf_M"] <- contador_medium
      Stats_Fusions[k, "Fusiones_conf_L"] <- contador_low
    } else {
      print(sprintf("Se encontraron 0 fusiones en %s", i))
    }


    Grupo <- Metadata[which(Metadata$ID == i), group]
    Stats_Fusions[k, "Grupo"] <- Grupo

    k = k+1

  }

  write.xlsx(Todos_FusionReport, file = sprintf("%s/Todos-FusionReports.xlsx", patients_dir))
  write.xlsx(Stats_Fusions, file = sprintf("%s/StatsFusions.xlsx", patients_dir))

  return(Todos_FusionReport, Stats_Fusions)

}
