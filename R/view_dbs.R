#' @title View Data Bases
#' @description Download the databases in excel file.
#' @export
#' @import openxlsx

view_dbs <- function() {
  #If a database was not created insert an empty data frame on the workbook
  path_mutdb <- sprintf("%s/mitorDB/DB/RDS_DB.rds", Sys.getenv('R_LIBS_USER'))

  if(file.exists(path_mutdb)) {
    mutDB <- readRDS(path_mutdb)
    mut_DB <- mutDB[[1]]
  } else {
    mut_DB <- data.frame()
  }

  path_cnvdb <- sprintf("%s/mitorDB/DB/cnvDB.RDS", Sys.getenv('R_LIBS_USER'))
  if(file.exists(path_cnvdb)) {
    cnv_DBs <- readRDS(path_cnvdb)
    ControlDB <- cnv_DBs[[1]]
    PatientsDB <- cnv_DBs[[2]]
    GenesDB <- cnv_DBs[[3]]
  } else {
    ControlDB <- data.frame()
    PatientsDB <- data.frame()
    GenesDB <- data.frame()
  }

  wb_dbs <- openxlsx::createWorkbook()

  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")

  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  addWorksheet(wb_dbs, "MutationsDB")
  addWorksheet(wb_dbs, "GenesDB")
  addWorksheet(wb_dbs, "PatientsDB")
  addWorksheet(wb_dbs, "ControlDB")

  writeData(wb_dbs, sheet = 1, x = mut_DB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")

  for (i in 1:nrow(mut_DB)){
    x <- mut_DB$Franklin[i]
    names(x) <- c("View on Franklin DB")
    class(x) <- "hyperlink"
    y <- mut_DB$VarSome[i]
    names(y) <- c("View on VarSome DB")
    class(y) <- "hyperlink"
    z <- mut_DB$dbSNP[i]
    names(z) <- c("View on dbSNP DB")
    class(z) <- "hyperlink"

    openxlsx::writeData(wb_dbs, sheet = 1, x = x, startRow = i+1, startCol = 13)
    openxlsx::writeData(wb_dbs, sheet = 1, x = y, startRow = i+1, startCol = 14)
    openxlsx::writeData(wb_dbs, sheet = 1, x = z, startRow = i+1, startCol = 9)
  }
  setColWidths(wb_dbs, sheet = 1, cols = c(1:ncol(mut_DB)), widths = "auto")

  writeData(wb_dbs, sheet = 2, x = GenesDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb_dbs, sheet = 3, x = PatientsDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")
  writeData(wb_dbs, sheet = 4, x = ControlDB,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")


  saveWorkbook(wb_dbs, paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DB/DataBases.xlsx", sep=""), overwrite = TRUE)
  browseURL(paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/DB/DataBases.xlsx", sep=""))

}

