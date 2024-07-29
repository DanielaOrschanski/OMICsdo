#' @title Update the MitoR-DataBase stats taken from HMTVAR
#' @description MitoR-DataBase takes information for each mutation from the HMTVAR DataBase. This DataBase is constantly updated.
#' We recommend to update the information of MitoR-DataBase at least once a month. With this function, the MitoR DataBase will be updated with the last
#' available information about each mutation from HMTVAR.
#' @param AddToPatients In case you want to update the information of the individual XLSX files from the patients as well, keep AddToPatients = TRUE.
#' Please be aware that this will take a while if you have many patients to update.
#' In case you want to update only one patient or just some of them, we recommend you to create the XLSX file from the VCF again. Check generateXLSX() for more information.
#' @return The DataBase XLSX file will automatically be updated. The same will happen with the patients XLSX files in case of AddToPatients = TRUE
#' @export
#' @examples
#' updateMitoRStats_HMTVAR()
#' updateMitoRStats_HMTVAR(addToPatients = TRUE)

updateMitoRStats_HMTVAR <- function(AddToPatients = FALSE){
  # Updates the MitoR DB from the HMTVAR database
  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))

  Freq_MitoR <- RDS_DB$Freq_MitoR[[1]]

  # Looks for the data given in the Freq_MitoR dataframe on the HMTVAR API
  updates <- buscarHMTVAR(Freq_MitoR)
  Freq_MitoR$clinVAR <- updates[[1]]
  Freq_MitoR$MitoMap <- updates[[2]]
  Freq_MitoR$dbSNP <- updates[[3]]
  Freq_MitoR$Omim <- updates[[4]]
  Freq_MitoR$Disease <- updates[[5]]
  Freq_MitoR$Freq_HMTVAR <- updates[[6]]

  # If it is not the first update, it deletes the last update row first
  updateRow <- data.frame(Paciente = "HMTVAR Update", Date = format(Sys.Date(), "%d/%m/%Y"), version_BWA = "-", version_GATK = "-", version_PICARD = "-", Last_Update = "-")
  if ("HMTVAR Update" %in% SoftData$Paciente) {
    SoftData <- SoftData[!SoftData$Paciente == "HMTVAR Update", ]
  }

  # Adds the update row
  SoftData <<- rbind(SoftData, updateRow)
  openxlsx::write.xlsx(mitoRStats, MitoRStats_Path, sheetName = "MitoRStats" , rowNames = FALSE)
  openxlsx::write.xlsx(SoftData, MitoRStats_Path, sheetName = "SoftData", rowNames = FALSE, overwrite = TRUE)

  if (AddToPatients == TRUE) {
    for (i in names(RDS_DB)[2:length(RDS_DB)]) {
      tryCatch(
        expr = {
          update_freq_MitoR(i)
        },
        error = function(e) {
          message(sprintf("An error occured while updating the patient number '%s'.
                          Please, try to update it individually by using the update_freq_MitoR() function.", i))

          print(e)
        },
        warning = function(w) {
          message(sprintf("An error occured while updating the patient number '%s'.
                          Please, try to update it individually by using the update_freq_MitoR() function.", i))
          print(w)
        },
        finally = {
          message(sprintf("Patient number '%s' has been succesfully updated.", i))
        }
      )
    }
  }
  RDS_DB$Freq_MitoR <- Freq_MitoR
  saveRDS(RDS_DB, sprintf("%s/RDS_DB.rds", mitor_db))
}


#' @title View Patient's Report
#' @description Downloads the xlsx report of a specific patient with the data updated.
#' @return message that indicates that the report has been updated and downloaded.
#' @export
#' @examples
#' view_patient()

update_freq_MitoR <- function(path_dir) {

  # Checks if the MitoR directory exists. If not the function is cancelled
  if (dir.exists(paste(path_dir, "/MitoR", sep = ""))) {
    file_list <- list.files(paste(path_dir, "/MitoR", sep=""))
    path_report <- file_list[endsWith(file_list, "MitoR_Report.xlsx")]
    patient_XLSX_file <- sprintf("%s/MitoR/%s", path_dir, path_report)
  } else {
    stop("This patient hasn't been analyzed. Please use Run_MitoR or select another patient.")
  }

  numPaciente <- basename(path_dir)

  # Adds or updates the MitoR frequencies of the patient analysis
  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))
  Freq_MitoR <- RDS_DB$Freq_MitoR

  # Reads the xlsx patient report
  paciente_SNP <- openxlsx::read.xlsx(patient_XLSX_file, sheet = 1)

  # Takes the mutations of the patient
  mutations <- RDS_DB[[sprintf("%s", numPaciente)]]["Mutations"]
  toAddFreq <- strsplit(unlist(mutations, use.names = FALSE), "/")

  # Save in a new variable the MitoR DB frequencies of the mutations found in the patient
  freqMut_MitoRStats <- c()
  for (i in (1:length(toAddFreq))){
    freqMut_MitoRStats <- c(freqMut_MitoRStats, Freq_MitoR[(Freq_MitoR$POS == as.integer(toAddFreq[[i]][2])) & (Freq_MitoR$ALT == toAddFreq[[i]][3]), ]$Freq_MitoR)
  }

  if ("Freq_MitoR" %in% colnames(paciente_SNP)) { # If its an update of the frequency
    paciente_SNP$Freq_MitoR <- freqMut_MitoRStats
  } else { # If this is the first time it adds the frequency
    paciente_SNP <- cbind(paciente_SNP, Freq_MitoR = freqMut_MitoRStats)
  }

  paciente_SNP_updated <- openxlsx::loadWorkbook(patient_XLSX_file)

  writeData(paciente_SNP_updated, sheet = "SNP-INDEL", x = paciente_SNP)

  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")
  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  writeData(paciente_SNP_updated, sheet = "SNP-INDEL", x = paciente_SNP,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")

  setColWidths(paciente_SNP_updated, sheet = "SNP-INDEL", cols = c(1:10, 12, 14, 15, 17:19), widths = "auto")

  for (i in 1:nrow(paciente_SNP)){
    x <- paciente_SNP$Franklin[i]
    names(x) <- c("View on Franklin DB")
    class(x) <- "hyperlink"
    y <- paciente_SNP$VarSome[i]
    names(y) <- c("View on VarSome DB")
    class(y) <- "hyperlink"
    z <- paciente_SNP$dbSNP[i]
    names(z) <- c("View on dbSNP DB")
    class(z) <- "hyperlink"

    openxlsx::writeData(paciente_SNP_updated, sheet = "SNP-INDEL", x = x, startRow = i+1, startCol = 14)
    openxlsx::writeData(paciente_SNP_updated, sheet = "SNP-INDEL", x = y, startRow = i+1, startCol = 15)
    openxlsx::writeData(paciente_SNP_updated, sheet = "SNP-INDEL", x = z, startRow = i+1, startCol = 10)
  }

  # Guardamos el actualizado en el mismo path que salio del selectDirectory
  saveWorkbook(paciente_SNP_updated, patient_XLSX_file, overwrite = TRUE)

  browseURL(patient_XLSX_file)
  return(sprintf("MitoR frequencies updated. Find the report at: %s", patient_XLSX_file))
}
