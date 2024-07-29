#' @title Convert VCF file to XLSX file
#' @description It generates a XLSX file with the information of the mutations available at the VCF file.
#' Creates a XLSX file from a VCF file. If you do not have one, we recommend you to see mitor::generateVCF()
#' By doing this, the genetist or the bioinformatitian is able to have the information on a more common file type.
#' Information about the mutation is also taken from HMTVAR. To use this function you must be online.
#'
#' The function has also the possibility to store the data in a database automatically generated. The database will also be a XLSX file.
#'
#' @param addToDataBase=TRUE While creating the XLSX file, the information is added to a DataBase.
#' If you want to keep add it as part of the analysis, it will be saved in the upper directory of the R1 and R2 files of the analyzed patient.
#' @return XLSX file
#' @import rjson
#' @export
#' @examples
#' generate_XLSX("../MitoR/patient1.vcf", addToDataBase=FALSE)
#' generate_XLSX("../MitoR/patient1.vcf")

generate_XLSX <- function(VCF_file, ref) {
  #mitor_sof <- sprintf("%s/mitorDB/Softwares", Sys.getenv('R_LIBS_USER'))
  #SAMTOOLS <- sprintf("%s/MitoRSoftware/Samtools/samtools-1.16.1/samtools", mitor_sof)
  #BWA <- sprintf("%s/MitoRSoftware/BWA/usr/bin/bwa", mitor_sof)
  #PICARD <- sprintf("%s/MitoRSoftware/picard-2.27.5/picard.jar", mitor_sof)
  #GATK <- sprintf("%s/MitoRSoftware/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
  #reference <- sprintf("%s/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_sof)

  #Agregado por mi:
  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

  #Checks all the softwares needed are downloaded ----
  linea_software <- grep("(?i)Samtools", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    Samtools <- downloadSamtools()
  } else {
    Samtools <-  sub("^Samtools ", "", linea_software)
  }

  linea_software <- grep("(?i)BWA", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    BWA <- downloadBWA()
  } else {
    BWA <-  sub("^BWA ", "", linea_software)
  }

  linea_software <- grep("(?i)PICARD", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    PICARD <- downloadPICARD()
  } else {
    PICARD <-  sub("^PICARD ", "", linea_software)
  }

  linea_software <- grep("(?i)GATK", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    GATK <- downloadGATK()
  } else {
    GATK <-  sub("^GATK ", "", linea_software)
  }

  if (ref == "HG38") {
    linea_software <- grep("(?i)HG38FASTA", softwares, ignore.case = TRUE, value = TRUE)
    if(length(nchar(linea_software)) == 0) {
      reference <- downloadHG38()
    } else {
      reference <- sub("^HG38FASTA ", "", linea_software)
    }
  } else {
    linea_software <- grep("(?i)HG19FASTA", softwares, ignore.case = TRUE, value = TRUE)
    if(length(nchar(linea_software)) == 0) {
      reference <- downloadHG19()
    } else {
      reference <- sub("^HG19FASTA ", "", linea_software)
    }
  }

  ##########################################################################
  numPaciente <- unlist(stringr::str_split(basename(VCF_file), "_"))[1]

  # Save the filtering parameters and take the VCF header out
  data_vcf <- analyzeVCF(VCF_file)

  # The filtering parameters will be added on the second sheet/WB
  params <- as.character(data_vcf[2])
  patient_vcf <- unlist(data_vcf[1])

  # Sets the colomn names: "CHROM" "POS" "ID" "REF" "ALT" "QUAL" "FILTER" "INFO" "FORMAT" "bar"
  col_names <- c(unlist(strsplit(patient_vcf[1], "\t")))
  col_names[1] <- "CHROM"

  # Sorts the VCF data: [1] "MT" [2] "73" [3] "." [4] "A" ...
  datos <- unlist(strsplit(patient_vcf[2:length(patient_vcf)], "\t"))

  # Merges the name of the colomns and the data into one matrix. Then creates a DF
  matriz <- t(matrix(datos, nrow = length(col_names), ncol = (length(patient_vcf) - 1)))
  filtered_DF <- as.data.frame(matriz, row.names = FALSE)
  colnames(filtered_DF) <- col_names


  # Saves the version of the used softwares
  versionGATK <- unlist(stringr::str_split(GATK, "[-/]"))[length(unlist(stringr::str_split(GATK, "[-/]")))-1]
  versionPICARD <- unlist(stringr::str_split(PICARD, "[/-]"))[length(unlist(stringr::str_split(PICARD, "[/-]")))-1]
  versionBWA <- "0.7.17"

  # Splits the information from the colomn bar: GT AD DP GQ PL
  #mio
  library(tidyr)
  vcf_df <- filtered_DF %>%
    separate(bar, into = c("GT", "AD", "DP", "GQ", "PL"), sep = ":")
  str(vcf_df)
  vcf_df$QUAL <- as.numeric(vcf_df$QUAL)
  #vcf_df <- vcf_df[which(vcf_df$QUAL > 100),]

  # From 0/0, 0/1, 1/1 to Hom, Het REF, Het ALT
  # This should only be applyed in nuclear genome. Not mitochondrial
  vcf_df$GT <- homHet(vcf_df$GT)
  Franklin <- sprintf("https://franklin.genoox.com/clinical-db/variant/snp/chr%s-%s-%s-%s", vcf_df$CHROM, vcf_df$POS, vcf_df$REF , vcf_df$ALT)
  VarSome <- sprintf("https://varsome.com/variant/hg38/%s%s%s%s%s%s%s?annotation-mode=germline", vcf_df$CHROM, "%3A", vcf_df$POS, "%3A", vcf_df$REF , "%3A", vcf_df$ALT)

  # Searching from the HMTVAR API information related to each mutation
  HMTVAR <- searchHMTVAR(vcf_df)
  clinvar <- HMTVAR[[1]]
  dbSNP <- sprintf("https://www.ncbi.nlm.nih.gov/snp/%s", HMTVAR[[2]])
  MitoMap <- HMTVAR[[3]]
  Omim <- HMTVAR[[4]]
  Disease <- HMTVAR[[5]]
  All_freq_h <- substr(HMTVAR[[6]], start = 0, stop = 5)

  # Delets unnecessary colomns
  #drop <- names(filtered_DF) %in% c("ID", "FILTER", "FORMAT", "bar")
  #filtered_DF <- filtered_DF[, !drop]

  # Merges the information from the HMTVAR API
  vcf_df$clinVAR <- clinvar
  vcf_df$dbSNP = dbSNP
  vcf_df$MitoMap = MitoMap
  vcf_df$Omim = Omim
  vcf_df$Disease = Disease
  vcf_df$Franklin = Franklin
  vcf_df$VarSome = VarSome
  vcf_df$Freq_HMTVAR = All_freq_h
  #vcf_df <- vcf_df[c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Allele_Depth",	"Deep_Coverage", "clinVAR",	"dbSNP",	"MitoMap",	"Omim",	"Disease", "Franklin", "VarSome", "INFO",	"GT",	"GQ",	"PL")]

  # Creates the DF used on the second sheet. Contains information about the patient ID, date of analysis,
  # filtering parameters and version of softwares
  filtered_DF2 <- data.frame(Patient = numPaciente, Date = format(Sys.Date(), "%d/%m/%Y"), Filter_Params = params,
                             VersionBWA = versionBWA, VersionGATK = versionGATK,
                             VersionPICARD = versionPICARD, Last_Update = format(Sys.Date(), "%d/%m/%Y"))

  # Path to the analysis XLSX file
  #XLSX_file <- paste(substr(VCF_file, start = 0, stop = (nchar(VCF_file)-4)), ".xlsx", sep = "", collapse = "")
  XLSX_file <- sprintf("%s/%s_OMICsdoVCF.xlsx", dirname(VCF_file), numPaciente)

  # Adds the patients information to the DB
    #ESTO ES LO QUE HAY QUE CHEQUEAR!!!!
  #freq_MitoR <- unlist(add_to_RDS(vcf_df, filtered_DF2, numPaciente))

  # Adds the information given by the MitoR DB
  #filtered_DF <- cbind(filtered_DF, Freq_MitoR = freq_MitoR)
  #filtered_DF <- filtered_DF[c("CHROM", "POS", "REF", "ALT", "Freq_MitoR", "Freq_HMTVAR",	"Allele_Depth",	"Deep_Coverage", "clinVAR",	"dbSNP",	"MitoMap",	"Omim",	"Disease", "Franklin", "VarSome", "INFO",	"GT",	"GQ",	"PL")]

  # Creates the WB with the patient's report DF and the software's report DF
  xlsx_SNP_Indels_report <- xlsx_report_design(vcf_df, filtered_DF2)

  return(list(xlsx_SNP_Indels_report[[1]], xlsx_SNP_Indels_report[[2]]))
  #return(list(vcf_df, filtered_DF2))
}


#' @title buscarHMTVAR
#' @description Takes information from the HMTVAR DataBase about the mutations we are interested in.
#' The taken data for each mutation is: clinvar, dbSNP, MitoMap, Omim, Disease, All_freq_h.
#' In case you are offline, the report will show an "Offline" value in the cells related to HmtVar
#' In case you the HmtVar website is down, the report will show an "HmtVar Down" value in the cells related to HmtVar
#' @param dataFrame It must have at least a REF, POS and ALT colomn name.
#' @return list with: clinvar, dbSNP, MitoMap, Omim, Disease, All_freq_h

searchHMTVAR <- function(dataFrame){
  clinvar <- c()
  dbSNP <- c()
  MitoMap <- c()
  Omim <- c()
  Disease <- c()
  All_freq_h <- c()
  # Checks if it has internet connection to enter the HMTVAR API page
  if (!httr::http_error(httr::GET("https://www.google.com"))){
    for (i in (1:nrow(dataFrame))){
      REF <- dataFrame[i, "REF"]
      POS <- dataFrame[i, "POS"]
      ALT <- dataFrame[i, "ALT"]

      # Limit time to connect to the HmtVar website (in seconds)
      limit_time <- 2

      time_of_response <- system.time ({
        datos <- system2("curl", sprintf("-X GET --max-time %s 'https://www.hmtvar.uniba.it/api/main/mutation/%s%s%s' -H 'accept: application/json' --insecure", limit_time, REF, POS, ALT), stdout = TRUE, wait = TRUE, timeout = limit_time + 1)
      })
      response_time <- time_of_response[[3]]
      print(response_time)

      # Checks if time limit was exceeded
      if (response_time > limit_time) {
        # Time limit exceeded
        clinvar <- "HmtVar Down"
        dbSNP <- "HmtVar Down"
        MitoMap <- "HmtVar Down"
        Omim <- "HmtVar Down"
        Disease <- "HmtVar Down"
        All_freq_h <- "HmtVar Down"
        break

      } else {
        # Time limit NOT exceeded
        datos1 <- rjson::fromJSON(datos)

        if (is.null(datos1$CrossRef$clinvar)){
          clinvar <- c(clinvar, "-")
        }
        else{
          clinvar <- c(clinvar, sprintf("%s", datos1$CrossRef$clinvar))
        }
        if (is.null(datos1$CrossRef$dbSNP)){
          dbSNP <- c(dbSNP, "-")
        }
        else{
          dbSNP <- c(dbSNP, sprintf("%s", datos1$CrossRef$dbSNP))
        }
        if (is.null(datos1$CrossRef$mitomap_associated_disease)){
          MitoMap <- c(MitoMap, "-")
        }
        else{
          MitoMap <- c(MitoMap, sprintf("%s", datos1$CrossRef$mitomap_associated_disease))
        }
        if (is.null(datos1$CrossRef$omim)){
          Omim <- c(Omim, "-")
        }
        else{
          Omim <- c(Omim, sprintf("%s", datos1$CrossRef$omim))
        }
        if (is.null(datos1$CrossRef$pubs_disease)){
          Disease <- c(Disease, "-")
        }
        else{
          Disease <- c(Disease, sprintf("%s", datos1$CrossRef$pubs_disease))
        }
        if (is.null(datos1$Variab$all_freq_h)){
          All_freq_h <- c(All_freq_h, "-")
        }
        else{
          All_freq_h <- c(All_freq_h, sprintf("%s", datos1$Variab$all_freq_h))
        }
      }
    }
  } else{
    clinvar <- "Offline"
    dbSNP <- "Offline"
    MitoMap <- "Offline"
    Omim <- "Offline"
    Disease <- "Offline"
    All_freq_h <- "Offline"
  }
  return(list(clinvar = clinvar, dbSNP = dbSNP, MitoMap = MitoMap, Omim = Omim, Disease = Disease, All_freq_h = All_freq_h))
}

#' @title homHet
#' @description Given the GATK parameters as 0/0, 0/1, 1/1, this function writes its meaning regarding the type of mutation.
#' @param GT List of GATK parameters as 0/0, 0/1, 1/1
#' @return The same list (if it is coming from the generateXLSX function, it is a dataFrame colomn) with the replaced meaning.
#' @examples homHet(list(0/0,0/1, 1/1))
homHet <- function(GT){
  for (i in (1:length(GT))){
    if (GT[i] == "0/0"){
      GT[i] <- "Hom Ref"
    } else if (GT[i] == "0/1"){
      GT[i] <- "Het"
    } else {
      GT[i] <- "Hom Alt"
    }
  }
  return(GT)
}

#' @title Modify the VCF file and delete the header
#' @description Discard the header of the VCF file and keep the variants information and the filter parameters (INDEL and SNP) used to generate the VCF
#' @param VCF_file Original VCF file
#' @return List of variant information, and filter parameters (INDEL and SNP)
#' @examples homHet(list(0/0,0/1, 1/1))

analyzeVCF <- function(VCF_file) {
  patient_vcf <- as.list(readLines(VCF_file))

  filter_params <- c(patient_vcf[grep("ID=mitor_indel_filter", patient_vcf, ignore.case = TRUE)], patient_vcf[grep("ID=mitor_snp_filter", patient_vcf)])
  filter_params <- paste("Indel Params = ", stringr::str_extract(filter_params[1], '(?<=Description=\\").*(?=\\")'), "SNP Params = ", stringr::str_extract(filter_params[2], '(?<=Description=\\").*(?=\\")'), collapse = " ")

  patient_vcf <- patient_vcf[-1:-(grep("#CHROM",patient_vcf) - 1)]
  return(list(patient_vcf, filter_params))
}



#' @title Design the MitoR report workbook
#' @description Designs the WB of SNP, INDEL and software reports. It creates the hyperlinks to Franklin, Varsome and dbSNP DB
#' @param filtered_DF DataFrame containing the SNP & INDEL report
#' @param filtered_DF2 DataFrame containing the software reports

xlsx_report_design <- function(filtered_DF, filtered_DF2) {

  wb_SNP_Indels <- openxlsx::createWorkbook()
  wb_Soft <- openxlsx::createWorkbook()

  COLNAMES_STYLE <- openxlsx::createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")

  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  openxlsx::addWorksheet(wb_SNP_Indels, "SNP_Indels")
  openxlsx::addWorksheet(wb_Soft, "Softwares")

  openxlsx::writeData(wb_SNP_Indels, sheet = 1, x = filtered_DF,
                      headerStyle = COLNAMES_STYLE,
                      borderStyle = "dashed",
                      borders = "columns", borderColour = "black")

  openxlsx::writeData(wb_Soft, sheet = 1, x = filtered_DF2,
                      headerStyle = COLNAMES_STYLE,
                      borderStyle = "dashed",
                      borders = "columns", borderColour = "black")

  openxlsx::setColWidths(wb_SNP_Indels, sheet = 1, cols = c(1:10, 12, 14, 15, 17:19), widths = "auto")
  openxlsx::setColWidths(wb_Soft, sheet = 1, cols = 1:ncol(filtered_DF2), widths = "auto")

  #hioervinculos -----------------------------------
  #for (i in 1:nrow(filtered_DF)){
  #  x <- filtered_DF$Franklin[i]
  #  names(x) <- c("View on Franklin DB")
  #  class(x) <- "hyperlink"
  #  y <- filtered_DF$VarSome[i]
  #  names(y) <- c("View on VarSome DB")
  #  class(y) <- "hyperlink"
  #  z <- filtered_DF$dbSNP[i]
  #  names(z) <- c("View on dbSNP DB")
  #  class(z) <- "hyperlink"

  #  openxlsx::writeData(wb_SNP_Indels, sheet = 1, x = x, startRow = i+1, startCol = 14)
  #  openxlsx::writeData(wb_SNP_Indels, sheet = 1, x = y, startRow = i+1, startCol = 15)
  #  openxlsx::writeData(wb_SNP_Indels, sheet = 1, x = z, startRow = i+1, startCol = 10)
  #}

  return(list(wb_SNP_Indels, wb_Soft))
}

