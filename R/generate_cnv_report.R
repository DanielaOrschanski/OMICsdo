#' @import openxlsx
#' @title Generation of report
#' @description Generates the workbook with the results of the patient's analysis.
#' You will find each CNV detected with it's information such as size, start, end,
#' association with Franklin and Varsome, etc.
#' @param CNV_calls dataframe generated using analyze_cnv function. It indicates the cnvs detected.
#' @param graph ggplot that represents the CNVs.
#' @param path_dir path of the directory of the patient that has been analized. It must contain R1 and R2.
#' @param transit transition probability: Transition probability of the hidden Markov Chain
#' from the normal copy number state to either a deletion or a duplication. Is used
#' in the generation of the CNV_calls indicated in the first parameter.
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' Is used in the generation of the CNV_calls indicated in the first parameter
#' @return wb workbook
#' @export

generate_cnv_report <- function(CNV_calls, graph = NULL, path_dir, transit, minoverlap, extra_bed = NULL) {

  CNV_calls <- generate_links(CNV_calls)

  id <- basename(path_dir)

  wb <- openxlsx::createWorkbook()

  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")


  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  addWorksheet(wb, paste(id, "CNVAnalisis"))
  #addWorksheet(wb, paste(id, "CNVplots"))
  addWorksheet(wb, paste(id, "Info"))
  writeData(wb, sheet = paste(id, "CNVAnalisis"), x = CNV_calls,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")


  setColWidths(wb, sheet = 1, cols = 1:4, widths = "auto")
  setColWidths(wb, sheet = 1, cols = 6:ncol(CNV_calls), widths = "auto")
  setColWidths(wb, sheet = 1, cols = 5, widths = 16)

  for(i in 1:nrow(CNV_calls)){
    x <- CNV_calls$Franklin[i]
    y <- CNV_calls$Varsome[i]
    names(x) <- c("View on Franklin DB")
    names(y) <- c("View on Varsome DB")
    class(x) <- "hyperlink"
    class(y) <- "hyperlink"

    indiceF <- which(colnames(CNV_calls) == "Franklin")
    indiceV <- which(colnames(CNV_calls) == "Varsome")

    writeData(wb, sheet = 1, x = x, startRow= i+1, startCol = indiceF)
    writeData(wb, sheet = 1, x = y, startRow= i+1, startCol = indiceV)

  }

  #Colour cells if the size of the variation is >500

  for (i in 1:nrow(CNV_calls)) {
    if(is.null(extra_bed)) {
      if (CNV_calls$BF[i] > 10) {
        addStyle(wb, 1, style = openxlsx::createStyle(fgFill = "thistle"), rows = i+1, cols = 1:ncol(CNV_calls))
      }
    } else {
      #tengo panel extra
      exons_list <- unlist(strsplit(CNV_calls$Exons[i], split=","))
      exons_list <- gsub(" ", "", exons_list)
      if (any(exons_list %in% extra_bed$name)) {
        addStyle(wb, 1, style = openxlsx::createStyle(fgFill = "thistle"), rows = i+1, cols = 1:ncol(CNV_calls))
      }
    }
  }

  #Register the parameters used to make the analysis:
  info1 <- paste("Minimum Overlap:", minoverlap, "  ||   Transition Probability:", transit, sep = " ")
  infodf <- data.frame("Patient" = id, "Date" = Sys.Date(), "VersionBWA"= "0.7.17",	"VersionGATK" = "4.3.0.0", "VersionPICARD" = "2.27.5", "CNV parameters" = info1)

  writeData(wb, sheet = paste(id, "Info"), x = infodf,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black",
            startRow = 1)

  setColWidths(wb, sheet = paste(id, "Info"), cols = 3:6, widths = "auto")

  #Add Plot
  #print(graph)
  #insertPlot(wb, 2, xy = c(1, 1), width = 25, height = 18, fileType = "png", units = "cm")

  #png("graph.png", width=1024, height=768, units="px", res=144)
  #print(graph)
  #dev.off()
  #insertImage(wb, 2, "graph.png", width=1024, height=768, units="px")



  return(wb)

}



#' @title Generation of Links
#' @description Generates hyperlinks that redirects the user to databases that provide more information about each variation, such as the pathogenicity and clinical implications.
#' @param CNV_calls dataframe with CNVs detected.
#' @return CNV_calls with 2 new columns "Franklin" and "Varsome".

generate_links <- function(CNV_calls) {
  #Order to generate links
  CNV_calls <- subset(CNV_calls, select = c("chromosome", "type", "size", "id", "ngenes", "nexons", "BF", "reads.expected", "reads.observed", "reads.ratio", "Genes","Exons"))
  colnames(CNV_calls)[1] <-"chr"

  for (i in 1:nrow(CNV_calls)) {
    mutF <- CNV_calls$id[i]
    mutF <- gsub(':', '-', mutF)
    mutV <- gsub("-", "%3A", mutF)
    addressF <- paste("https://franklin.genoox.com/clinical-db/variant/sv/", mutF, sep = "")
    addressV <- paste("https://varsome.com/cnv/hg19/", mutV, sep = "")

    if (CNV_calls$type[i] == "deletion") {
      addressF <- paste(addressF, "-DEL", sep = "")
      addressV <- paste(addressV, "%3ADEL?", sep = "")
    } else {
      addressF <- paste(addressF, "-DUP", sep = "")
      addressV <- paste(addressV, "%3ADUP?", sep = "")
    }

    CNV_calls$Franklin[i] <- addressF
    CNV_calls$Varsome[i] <- addressV

  }

  #CNV_calls <- CNV_calls[, c(1,2,3,4,5,6,11,12,7,8,9,10)]
  return(CNV_calls)
}
