#' @import ExomeDepth
#' @import tibble
#' @import Rsamtools
#' @title Add patient to Patients DB
#' @description Saves reads' counts of each gene from the sequence of the patient in a dataset.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @return PatientsDB and Pcounts
#' @export
#' @examples PatientsDB <-  AddPatient( path_dir ="/home/sam/Patients/123/123.fastq")

add_patient_cnv <- function(path_dir, bed, ref) {

  initial <- initialize_db(bed = bed, ref = ref)
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]

  id <- basename(path_dir)

  if (id %in% colnames(PatientsDB)){ #If the patients is in Patients DB
    message("This patient is already in your database")
    return(PatientsDB)

  } else if (id %in% colnames(ControlDB)) {
    indic <- which(colnames(ControlDB) == id)
    Pcounts <- ControlDB[, indic]
    message("This patient is already in your control database")

  } else { #If the patient is not on any DB
    message("This patient is not already in the data base. This process will take a few minutes")

    #Check if there are BAM or BAI files already generated or generates them:
    bam.control <- getBAMBAI(path_dir, ref = ref)
    #bam.control <- "/media/16TBDisk/Daniela/CNVs/MuestrasCNVs/37513/OMICsdo/37513_recal.bam"

    include <- ifelse(ref == "HG38", TRUE, FALSE)
    cts <- getBamCounts(bed.frame = bed,
                        bam.files = bam.control ,
                        include.chr = include, #if set to TRUE, will add the string 'chr' to the chromosome names of the target BED file.
                        referenceFasta = NULL)

    Pcounts <- data.frame(Pcounts = cts[, ncol(cts)])

  }

  #Save counts on PatientsDB
  if (max(dim(PatientsDB)) == 0) { #this will be the first patient in DB
    PatientsDB <- data.frame(chromosome = bed[,1], start = bed[,2], end = bed[,3], name=bed[,4] , Pcounts = cts[, ncol(cts)])
    names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
    #colnames(PatientsDB)[ncol(PatientsDB)] <- id
    message("Your first patient has been saved!")

  } else { #not the first control patient
    PatientsDB <- add_column(PatientsDB, Pcounts)
    names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
    #colnames(PatientsDB)[ncol(PatientsDB)] <- id
    message("Another patient has been saved!")
  }

  #Save the updated Patients DB
  #libPath <- dirname(dirname(dirname(system.file(package = "MitoR"))))
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/OMICsdo/cnvDB.RDS", sep=""))
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))

  #return(list(PatientsDB, Pcounts))
  return(PatientsDB)
}
