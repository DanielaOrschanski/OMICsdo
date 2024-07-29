#' @import org.Hs.eg.db
#' @title Analyze copy number variations (CNVs) of a patient sample.
#' @description Gjrdbhjdbj
#' @param path_fc path to the feature counts rds file
#' @param path_metadata path to the metadata excel file
#' @return a list:
#' @export
#library(org.Hs.eg.db)

FCtoDGEList <- function(path_fc, path_metadata) {
  #path_fc <- "/media/4tb1/Daniela/Environ/Fibroblastos/MuestrasFibroblastos/Colon/FeatureCount_Report-Fibroblastos.rds"
  #path_metadata <- "/media/4tb1/Daniela/Environ/RedNeuronal/Metadata Environ -  Mama.xlsx"

  FeatureCount_Report <- readRDS(path_fc)
  counts <- FeatureCount_Report$counts
  nombres <- strsplit(colnames(counts), split="_")
  colnames(counts) <- sapply(nombres, function(x) x[1])

  geneid <- rownames(counts)
  genes <- mapIds(org.Hs.eg.db, geneid, "SYMBOL", "ENTREZID")
  environ_genes <- as.data.frame(cbind(geneid, genes))
  colnames(environ_genes) <- c("ENTREZID","SYMBOL")

  Metadata_Mama_Fibro <- read_excel(path_metadata)
  rownames(Metadata_Mama_Fibro) <- Metadata_Mama_Fibro$`ID ENVIRON`


  #which(colnames(counts) %in% Metadata_Mama_Fibro$`ID ENVIRON`)
  rownames(Metadata_Mama_Fibro) == colnames(counts)

  Br_Fibroblastos <- DGEList(counts,
                             lib.size = colSums(counts),
                             norm.factors = calcNormFactors(counts),
                             samples = Metadata_Mama_Fibro,
                             genes = environ_genes)

  rownames(counts) <- environ_genes[,2]

  #write.xlsx(counts, file = sprintf("%s/Counts.xlsx", dirname(path_fc)))
  write_rds(counts, file = sprintf("%s/Counts.rds", dirname(path_fc)))
  #write.xlsx(counts, file = sprintf("%s/DGEList.rds", dirname(path_fc)))
}

