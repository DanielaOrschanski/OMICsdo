#' @title Analyze mutations and variations
#' @description Generates the complete analysis to detect SNPs, indels and CNVs.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain
#' from the normal copy number state to either a deletion or a duplication.
#' @return a list: Patients DB (updated if the patient that has been analyzed is new),
#' CNV calls (dataframe with information about CNVs),
#' GenesDB (updated data base with information about the genes that has been detected
#' as variant)
#' @export
#' @examples
#' outAnalisis <- MitoRAnalysis(path_dir = "~/PacientesMito/46608/BM22-46608_R1.fastq.gz")
#' PatientsDB<-outAnalisis[[1]]
#' CNV_calls<-outAnalisis[[2]]
#' all.exons<-outAnalisis[[3]]

mitor_analysis <- function(path_dir, cnv_minoverlap = 0.0001, cnv_transit = 1, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {

  id <- basename(path_dir)

  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #         CNVs ANALYSIS
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  out_cnv <- detection_cnv(path_dir, cnv_minoverlap, cnv_transit)
  CNV_calls <- out_cnv[[1]]
  CNV_calls <- generate_links(CNV_calls)
  graph <- out_cnv[[2]]
  #wb_cnv <- generate_cnv_report(CNV_calls, graph, path_dir, transit = cnv_transit, minoverlap = cnv_minoverlap)

  message("CNVs' analysis has been done!")

  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  #         SNP-INDEL ANALYSIS
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  workbooks <- SNP_Indel_report(path_dir, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  message("SNPs&INDELs' analysis has been done!")

  #workbooks <- list(wb_SNP_Indel_patient, wb_SNP_Indel_soft, plots_to_xlsx)
  wb1<- workbooks[[1]]
  wb2<- workbooks[[2]]
  plots_to_xlsx <- workbooks[[3]]

  #file.remove(sprintf("%s/MitoR/%s_sortedR1R2.bam", path_dir, id))
  #file.remove(sprintf("%s/MitoR/%s_sortedR1R2.bam.bai", path_dir, id))

  # Create complete workbook
  wb_complete <- createWorkbook()

  #Format features:
  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")
  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  # Mutation's Sheet
  library(openxlsx)
  #d1 <- read.xlsx(wb1, sheet = 1)
  if(class(wb1)== "data.frame") {
    d1 <- wb1
  } else {
    d1 <- read.xlsx(wb1, sheet = 1)
  }

  addWorksheet(wb_complete, sheetName = "SNP-INDEL")
  writeData(wb_complete, sheet = "SNP-INDEL", x = d1)

  writeData(wb_complete, sheet = "SNP-INDEL", x = d1,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")

  setColWidths(wb_complete, sheet = "SNP-INDEL", cols = c(1:10, 12, 14, 15, 17:19), widths = "auto")

  for (i in 1:nrow(d1)){
    x <- d1$Franklin[i]
    names(x) <- c("View on Franklin DB")
    class(x) <- "hyperlink"
    y <- d1$VarSome[i]
    names(y) <- c("View on VarSome DB")
    class(y) <- "hyperlink"
    z <- d1$dbSNP[i]
    names(z) <- c("View on dbSNP DB")
    class(z) <- "hyperlink"

    openxlsx::writeData(wb_complete, sheet = "SNP-INDEL", x = x, startRow = i+1, startCol = 14)
    openxlsx::writeData(wb_complete, sheet = "SNP-INDEL", x = y, startRow = i+1, startCol = 15)
    openxlsx::writeData(wb_complete, sheet = "SNP-INDEL", x = z, startRow = i+1, startCol = 10)
  }

  #Mutation's Plot
  addWorksheet(wb_complete, sheetName = "Plot_Mut")
  for (i in 1:length(plots_to_xlsx)) { # i is related to the X axis
    for (j in 1:length(plots_to_xlsx[[i]])) { # j is related to the Y axis
      print(plots_to_xlsx[[i]][[j]])
      replayPlot(plots_to_xlsx[[i]][[j]])
      openxlsx::insertPlot(wb_complete, 2, xy = c((i*9)-8, (j*23)-21), width = 18.5, height = 12, fileType = "png", units = "cm")

      #png(sprintf("%s/MitoR/graph_mut.png", path_dir), width = 18.5, height = 12, units = "cm", res = 144)
      #print(plots_to_xlsx[[i]][j])
      #dev.off()

      #if (file.exists()) {
      #  insertImage(wb_complete, 2, sprintf("%s/MitoR/graph_mut.png", path_dir), width = 18.5, height = 12, startRow = (i*9)-8, startCol = (j*23)-21, units = "cm")
      #  system2("rm", sprintf("%s/MitoR/graph_mut.png", path_dir))
      #}
    }
  }

  #CNV's Sheet
  #CNV_calls <- read.xlsx(wb_cnv, sheet = 1)
  addWorksheet(wb_complete, sheetName = "CNV")
  writeData(wb_complete, sheet = "CNV", x = CNV_calls,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")

  setColWidths(wb_complete, sheet = "CNV", cols = 1:4, widths = "auto")
  setColWidths(wb_complete, sheet = "CNV", cols = 6:ncol(CNV_calls), widths = "auto")
  setColWidths(wb_complete, sheet = "CNV", cols = 5, widths = 16)

  for(i in 1:nrow(CNV_calls)){
    x <- CNV_calls$Franklin[i]
    y <- CNV_calls$Varsome[i]
    names(x) <- c("View on Franklin DB")
    names(y) <- c("View on Varsome DB")
    class(x) <- "hyperlink"
    class(y) <- "hyperlink"

    indiceF <- which(colnames(CNV_calls) == "Franklin")
    indiceV <- which(colnames(CNV_calls) == "Varsome")

    writeData(wb_complete, sheet = "CNV", x = x, startRow= i+1, startCol = indiceF)
    writeData(wb_complete, sheet = "CNV", x = y, startRow= i+1, startCol = indiceV)
  }

  #Colour cells if the size of the variation is >500
  for (i in 1:nrow(CNV_calls)) {
    if (!is.na(CNV_calls$size[i]) && CNV_calls$size[i] > 500) {
      addStyle(wb_complete, sheet = "CNV" , style = openxlsx::createStyle(fgFill = "thistle"), rows = i+1, cols = 1:ncol(CNV_calls))
    }
  }

  #CNV's Visualization Sheet
  addWorksheet(wb_complete, sheetName = "Plot_CNV")
  print(graph)
  insertPlot(wb_complete, "Plot_CNV", xy = c(1, 1), width = 25, height = 18, fileType = "png", units = "cm")
  #png(sprintf("%s/MitoR/graph_cnv.png", path_dir), width=1024, height=768, units="px", res=144)
  #print(graph)
  #dev.off()

  #if (file.exists()) {
  #  insertImage(wb_complete, 2, sprintf("%s/MitoR/graph_cnv.png", path_dir), width = 18.5, height = 12, startRow = (i*9)-8, startCol = (j*23)-21, units = "cm")
  #  system2("rm", sprintf("%s/MitoR/graph_cnv.png", path_dir))
  #}

  #Software's Information Sheet
  if(class(wb2) == "data.frame") {
    d3 <- wb2
  } else {
    d3 <- read.xlsx(wb2, sheet = 1)
  }

  addWorksheet(wb_complete, sheetName = "Soft")
  writeData(wb_complete, sheet = "Soft", x = d3,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns",
            borderColour = "black")

  info_cnv <- data.frame("CNV parameters" = paste("Minimum Overlap:", cnv_minoverlap, "  ||   Transition Probability:", cnv_transit, sep = " "))

  writeData(wb_complete, sheet = "Soft", x = info_cnv,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns",
            borderColour = "black",
            startCol = ncol(d3) + 1)

  setColWidths(wb_complete, sheet = "Soft", cols = 1:ncol(d3) + 1, widths = "auto")

  openxlsx::saveWorkbook(wb_complete, paste(path_dir,"/MitoR/", id, "_MitoR_Report.xlsx", sep = ""),  overwrite = TRUE)
  browseURL(paste(path_dir,"/MitoR/", id, "_MitoR_Report.xlsx", sep = ""))

}
