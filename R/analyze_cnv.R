#' @import ExomeDepth
#' @import tibble
#' @import dplyr
#' @title Analyze copy number variations (CNVs) of a patient sample.
#' @description Generate an excel report with 2 sheets: in the first one there is
#' a table with information about the CNVs detected such as location, size, and hyperlinks to online data bases.
#' On the second sheet there is a visualization of the CNVs within the whole mitochondrial genome of the patient.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain from the normal copy number state to either a deletion or a duplication.
#' This parameter is taken from ExomeDepth algorithms.
#' @return a list:
#' - Patients DB (updated if the patient that has been analyzed is new).
#' - CNV calls (dataframe with information about CNVs).
#' - GenesDB (updated data base with information about the genes that has been detected as variant.
#' - Graph with the duplications colored in green and deletions colored in red.
#' @export
#' @examples
#' outAnalisis <- analyze_cnv(path_dir = "~/Patients/46608", transit = 0.7)
#' PatientsDB <- outAnalisis[[1]]
#' CNV_calls <- outAnalisis[[2]]
#' GenesDB <- outAnalisis[[3]]
#' graph <- outAnalisis[[4]]

analyze_cnv <- function(path_dir, bed, ref, extra_bed = NULL, minoverlap = 0.5, transit = 1, Nspan = 0.8) {

  out_cnv <- detection_cnv(path_dir = path_dir, bed = bed, ref = ref, minoverlap = minoverlap, transit = transit)
  CNV_calls_f <- out_cnv[[1]]
  CNV_calls_filt <- out_cnv[[2]]
  all.exons <- out_cnv[[3]]
  #CNV_calls <- out_cnv[[1]]
  #graph <- out_cnv[[2]]

  #if(ref == "HG19" & is.null(bed) ) {
    #Usamos el bed completo por default
  #  data("exons.hg19")
  #  bed <<- exons.hg19
  #}


  id <- basename(path_dir)


  #wb <- generate_cnv_report(CNV_calls, graph = graph, path_dir, transit = transit, minoverlap = minoverlap)
  #CNV_calls_f <- na.omit(CNV_calls)
  #CNV_calls_filt <- CNV_calls_f[CNV_calls_f$size>1000000,]

  #wb <- generate_cnv_report(CNV_calls_filt, path_dir = path_dir, transit = transit, minoverlap = minoverlap)
  cnv_graph <- plot_cnv(all.exons = all.exons, patient_dir = path_dir, Nspan = Nspan)

  wb <- generate_cnv_report(CNV_calls_f, path_dir = path_dir, transit = transit, minoverlap = minoverlap, extra_bed= extra_bed)


  if (file.exists(paste(path_dir, "/OMICsdo", sep=""))) {
    openxlsx::saveWorkbook(wb, paste(path_dir, "/OMICsdo/", id, "_CNVreport.xlsx", sep = ""),  overwrite = TRUE)
    #browseURL(paste(path_dir, "/OMICsdo/", id, "_CNVreport.xlsx", sep = ""))
  } else {
    openxlsx::saveWorkbook(wb, paste(path_dir, "/", id, "_CNVreport.xlsx", sep = ""),  overwrite = TRUE)
    #browseURL(paste(path_dir, "/", id, "_CNVreport.xlsx", sep = ""))
  }

  message("Report has been succesfully generated")

  initial <- initialize_db()
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]

  return(list(CNV_calls_f, GenesDB, PatientsDB, all.exons, cnv_graph))
  #return(list(PatientsDB, CNV_calls, GenesDB, graph))
}


#' @title Detect copy number variations (CNVs) of a patient sample.
#' @description Generates a dataframe and a visual representation of the CNVs detected.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain from the normal copy number state to either a deletion or a duplication.
#' This parameter is taken from ExomeDepth algorithms.
#' @return a list:
#' - CNV calls: dataframe with information about CNVs.
#' - Graph: visual representation of CNVs. Duplications colored in green and deletions colored in red.
#' @import dplyr

detection_cnv <- function(bed, ref, path_dir, minoverlap, transit) {

  initial <- initialize_db(bed = bed, ref = ref)
  ControlDB <- initial[[1]]
  PatientsDB <- initial[[2]]
  GenesDB <- initial[[3]]

  #data("exons.hg19")
  #bed <<- exons.hg19

  if (nrow(ControlDB) == 0) {
    return("You must have at least one patient in Control DB. Use command add_control_cnv().")
  }

  id <- basename(path_dir)

  # Avoid using the patient as a reference for itself
  initialCDB <- ControlDB # save the previous controlDB complete

  if(id %in% colnames(ControlDB)) {
    indic <- which(colnames(ControlDB) == id)
    Pcounts <- ControlDB[, indic]

    #Delete the patient from ControlDB
    ControlDB <- subset(ControlDB, select = -indic)
    if (max(dim(ControlDB)) == 0 | ncol(ControlDB) == 0) { stop("You need to have at least one control patient") }

    #Add the patient to PatientsDB if its not there
    if (!(id %in% colnames(PatientsDB))) {
      if (max(dim(PatientsDB)) == 0) { #this will be the first patient in DB
        PatientsDB <- data.frame(Pcounts)
        names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
        message("Your first patient has been saved!")

      } else { #not the first patient on the DB
        PatientsDB <- add_column(PatientsDB, Pcounts)
        names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
        message("Another patient has been saved!")
      }
    }
  }

  if (!(id %in% colnames(PatientsDB))) { #If its a new patient---------------------------------------------------------
    message("This patients is not already on your database.The process may take a few minutes")

    #Check if there are BAM or BAI files already generated or generates them:
    bam.control <- getBAMBAI(path_dir, ref = ref)

    #Intenta con TRUE el include.chr y si da error que intente con FALSE
    tryCatch(
      {
        #include <- ifelse(ref == "HG38", FALSE, TRUE)
        include <- FALSE
        cts <<- getBamCounts(bed.frame = bed,
                            bam.files = bam.control ,
                            include.chr = include, #if set to TRUE, this function will add the string 'chr' to the chromosome names of the target BED file.
                            referenceFasta = NULL)
      },
      error = function(e) {
        cts <<- getBamCounts(bed.frame = bed,
                            bam.files = bam.control ,
                            include.chr = T, #if set to TRUE, this function will add the string 'chr' to the chromosome names of the target BED file.
                            referenceFasta = NULL)
      }
    )

    Pcounts <- data.frame(Pcounts = cts[, ncol(cts)])

    #Add the new patient to the PatientsDB
    if (max(dim(PatientsDB)) == 0) { #this will be the first patient in DB
      PatientsDB <- data.frame(Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Your first patient has been saved!")

    } else { #not the first control patient
      PatientsDB <- add_column(PatientsDB, Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Another patient has been saved!")
    }
  }#------------------------------------------------------------------------------------------------------------------

  tryCatch(
    {
      cnvs_out <- cnv_from_counts(id = id, bed = bed, ControlDB = ControlDB, PatientsDB = PatientsDB, minoverlap = minoverlap, transit = transit)
      CNV_calls <- cnvs_out[[1]]
      all.exons <- cnvs_out[[2]]
      print("There are CNVs detected")
    },
    error = function(e) {
      stop("No CNVs were found. Try again using a higher value of transit")
    }
  )

  #CNV_calls_f <- na.omit(CNV_calls)
  #CNV_calls_filt <- CNV_calls_f[CNV_calls_f$size > 10000,]

  #Me quedo con las primeras 20 CNVs (más grandes de tamaño) para graficarlas
  CNV_calls_filt <- CNV_calls[order(CNV_calls$size, decreasing = TRUE),]
  CNV_calls_filt <- head(CNV_calls_filt, 20)
  #graph <- plot_cnv(CNV_calls_filt)
  message("The graph has been succesfully generated")

  #Aca viene el problema:
  #GenesDB <- update_DB_cnv(CNV_calls = CNV_calls, id, ControlDB, PatientsDB, GenesDB)

  ControlDB <- initialCDB #if the patient was on ControlDB it has to be removed just for the analysis. ##
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/OMICsdoSof/cnvDB.RDS", sep=""))
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(dirname(system.file(package = "OMICsdo")), "/OMICsdoSof/cnvDB.RDS", sep=""))

  message("Genes DB has been succesfully updated")

  return(list(CNV_calls, CNV_calls_filt, all.exons))
  #return(list(CNV_calls, graph))
}


#' @title Add another patient's information to Genes DB
#' @description Adds the information about the genes affected by the CNVs of a new patient
#' and keeps track of the frequency and type of variation of each gene.
#' @param CNV_calls dataframe which contains the information about the CNVs of the new patient.
#' @param id number or character that indicates which patient is.
#' @param ControlDB database of control patients.
#' @param PatientsDB database of patients to analyze or analyzed.
#' @param GenesDB database of genes affected.
#' @return Genes DB updated.
#' @import dplyr
update_DB_cnv <- function(CNV_calls, id, ControlDB, PatientsDB, GenesDB) {

  #If the patient was analyzed before, delete the old data:
  if(id %in% colnames(GenesDB)) {
    indic <- which(colnames(GenesDB) == id)
    GenesDB <- subset(GenesDB, select = -indic)
  }

  #x <- c()
  #e_var <- c()
  #for (i in CNV_calls$name.gene){
  #  x <- append(x, strsplit(i, split = ","))
  #}
  e_var <- unlist(strsplit(CNV_calls$Genes, split = ","))

  #for (l in 1:length(x)){
  #  for (g in 1:length(x[[l]])){
  #    e_var <- append(e_var, x[[l]][g])
  #  }
  #}

  Complete <- data.frame("Gene" = e_var)

  #Amount of duplications or delections for each gene
  all_dup <- c()
  all_del <- c()

  #Create 2 lists: one for the duplicated genes and another with the deleted
  g_dup <- unlist(CNV_calls[CNV_calls$type == "duplication", "Genes" ])
  if (!(length(g_dup) == 0)) {
    all_dup <- c()
    q <- c()

    for (i in g_dup) {  q <- append(q, strsplit(i, split =","))  }
    for (l in 1:length(q)) {
      for (g in 1:length(q[[l]])) { all_dup <- append(all_dup, q[[l]][g]) }
    }
  }

  g_del <- unlist(CNV_calls[CNV_calls$type == "deletion", "name.gene" ])
  if (!(length(g_del) == 0)) {
    all_del <- c()
    w <- c()

    for (i in g_del) { w <- append(w, strsplit(i, split=",")) }
    for (l in 1:length(w)){
      for (g in 1:length(w[[l]])) { all_del<-append(all_del, w[[l]][g]) }
    }
  }

  #Indicate if the genes was duplicated or deleted for each patient (column):
  Complete$Paciente[Complete$Gene %in% all_dup ] <- "duplication"
  Complete$Paciente[Complete$Gene %in% all_del ] <- "deletion"


  if (ncol(GenesDB) == 1) { #If its the first patient in Genes DB
    GenesDB <- data.frame(Complete)
    GenesDB$Frec_dup <- 0
    GenesDB$Frec_del <- 0
    message("The first patient has been saved in GenesDB")

  } else {
    p <- as.data.frame(Complete[, c(1, 2)])
    GenesDB <- merge(GenesDB, p, by = "Gene", all = T)
    GenesDB <- select(GenesDB, -c("Frec_dup", "Frec_del", "Frec_variant") )
    message("Another patient has been saved in GenesDB")
  }

  names(GenesDB)[names(GenesDB) == "Paciente"] <- id

  GenesDB$Frec_dup <- 0
  GenesDB$Frec_del <- 0

  for (i in 1:nrow(GenesDB)) {
    for (j in 2:(ncol(GenesDB) - 2)) {
      if (!(is.na(GenesDB[i,j])) & GenesDB[i, j] == "deletion" ) {
        GenesDB$Frec_del[i] <- GenesDB$Frec_del[i] + 1
      } else if (!(is.na(GenesDB[i,j])) & GenesDB[i, j] == "duplication") {
        GenesDB$Frec_dup[i] <- GenesDB$Frec_dup[i] + 1
      }
    }
  }

  #Total amount of variations per gene:
  GenesDB$Frec_variant <- GenesDB$Frec_dup + GenesDB$Frec_del
  #libPath <- dirname(dirname(dirname(system.file(package = "OMICsdo"))))
  saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/OMICsdoSof/cnvDB.RDS", sep=""))

  return(GenesDB)

}
