#' @title Initialize Data Bases
#' @description Generation of 3 empty databases: Control(normal references),
#'  Patients (analyzed) and Genes (genes detected with variation) if it is the first
#'  time that the user uses this package. Else, it reads the databases already
#'  created.
#' @import ExomeDepth
#' @return list with 3 databases
#' @export
#' @examples
#' initial <- InitializeDB()
#' ControlDB <- initial[[1]]
#' PatientsDB <- initial[[2]]
#' GenesDB <- initial[[3]]

initialize_db <- function(bed = NULL, ref = "HG19") {
  #library(ExomeDepth)

  if(ref == "HG19" & is.null(bed)) {
    data("exons.hg19")
    bed <<- exons.hg19
  }

  bed <<- bed
  ref<<- ref


  #if (!(file.exists(sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))))) {
  if (!(file.exists(sprintf("%s/OMICsdoSof/cnvDB.RDS", Sys.getenv('R_LIBS_USER'))))) {

    wb1 <- data.frame()
    wb2 <- data.frame()
    wb3 <- data.frame("Gene" = bed[, 4])


    if(file.exists(paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))) {
      saveRDS(list(wb1, wb2, wb3), file = paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))
    } else if (file.exists(paste(Sys.getenv('R_LIBS_USER'), "/OMICsdoSof/cnvDB.RDS", sep=""))) {
      saveRDS(list(wb1, wb2, wb3), file = paste(Sys.getenv('R_LIBS_USER'), "/OMICsdoSof/cnvDB.RDS", sep=""))
    } else {
      #No existe el cnvDB.RDS
      saveRDS(list(wb1, wb2, wb3), file = paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))
    }

    message("Data Bases have been created. You can find them in OMICsdoSof folder")
  }

  #cnv_DBs <- readRDS(paste(Sys.getenv('R_LIBS_USER'), "/OMICsdoSof/cnvDB.RDS", sep=""))
  cnv_DBs <- readRDS(paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))
  ControlDB <- cnv_DBs[[1]]
  PatientsDB <- cnv_DBs[[2]]
  GenesDB <- cnv_DBs[[3]]

  return(list(ControlDB, PatientsDB, GenesDB, bed, ref))

}



