#' @importFrom ExomeDepth getBamCounts
#' @import tibble
#' @import Rsamtools
#' @title Add patient to Control DB
#' @description Saves reads' counts of each gene from the sequence of the patient in a reference data set use as a control data base.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta/q files).
#' @return Control DB updated
#' @examples ControlDB <- AddControl(path_dir ="/home/sam/Patients/123/123.fastq")
#' @export

add_control_cnv <- function(path_dir, bed, ref) {

  initial <- initialize_db(bed = bed, ref = ref)
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]

  id <- basename(path_dir)

  if (!(id %in% colnames(ControlDB))) { #Make sure the patient is not already in controlDB

    if (id %in% colnames(PatientsDB)) { #If the patients is in PatientsDB-----------------------------
      message("This patient is already in the data base. It will be copied to the control data base")
      indic <- which(colnames(PatientsDB) == id)
      Pcounts <- PatientsDB[, indic]

    } else { #If the patient is not on the DB------------------------------------------------------------
      message("This patient is not already in the data base. This process will take a few minutes")

      #Check if there are BAM or BAI files already generated or generates them:
      bam.control <- getBAMBAI(path_dir, ref = ref)

      include <- ifelse(ref == "HG38", TRUE, FALSE)
      cts <- getBamCounts(bed.frame = bed,
                          bam.files = bam.control ,
                          include.chr = include, # if set to TRUE, will add the string 'chr' to the chromosome names of the target BED file.
                          referenceFasta = NULL)

      Pcounts <- data.frame(Pcounts = cts[, ncol(cts)])

    }#---------------------------------------------------------------------------

    #update the database

    if (max(dim(ControlDB)) == 0) { #this will be the first control patient
      #ControlDB <- data.frame(Pcounts)
      ControlDB <- data.frame(chromosome = bed[,1], start = bed[,2], end = bed[,3], name=bed[,4] , Pcounts = cts[, ncol(cts)])
      names(ControlDB)[names(ControlDB) == "Pcounts"] <- id
      message("Your first control patient has been saved!")

    } else { #not the first control patient
      ControlDB <- add_column(ControlDB, Pcounts)
      names(ControlDB)[names(ControlDB) == "Pcounts"] <- id
      message("Another control patient has been saved!")
    }

    #libPath <- dirname(dirname(dirname(system.file(package = "MitoR"))))
    saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/OMICsdoSof/cnvDB.RDS", sep=""))
    saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))

  } else {
    message("This patient is already on your Control Data Base")
    return( ControlDB)
  }

  return(ControlDB)

}

