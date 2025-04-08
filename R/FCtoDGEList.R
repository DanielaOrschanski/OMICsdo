#' @import org.Hs.eg.db
#' @title Analyze copy number variations (CNVs) of a patient sample.
#' @description Gjrdbhjdbj
#' @param path_fc path to the feature counts rds file
#' @param path_metadata path to the metadata excel file
#' @return a list:
#' @export
#library(org.Hs.eg.db)

FCtoDGEList <- function(path_fc, path_metadata) {
  #path_fc <- "/media/4tb2/Daniela/MuestrasChile2/Muestras/FeatureCount_Report.rds"
  #path_metadata <- "/media/4tb2/Daniela/MuestrasChile2/Meta_data_TPJ_RNAseqJuly31- ACTUALIZADA.xlsx"

  FeatureCount_Report <- readRDS(path_fc)
  counts <- FeatureCount_Report$counts
  nombres <- strsplit(colnames(counts), split="_Aligned")
  colnames(counts) <- sapply(nombres, function(x) x[1])

  geneid <- rownames(counts)
  genes <- mapIds(org.Hs.eg.db, geneid, "SYMBOL", "ENTREZID")
  environ_genes <- as.data.frame(cbind(geneid, genes))
  colnames(environ_genes) <- c("ENTREZID","SYMBOL")

  Metadata_Mama_Fibro <- as.data.frame(read_excel(path_metadata))
  columnas_nombres_ids <- which(sapply(Metadata_Mama_Fibro, function(x) any(colnames(counts) %in% x)))[1]
  colnames(Metadata_Mama_Fibro)[columnas_nombres_ids] <- "ID"
  Metadata_Mama_Fibro <- Metadata_Mama_Fibro[which((Metadata_Mama_Fibro$ID %in% colnames(counts))),]

  Metadata_Mama_Fibro <- Metadata_Mama_Fibro[order(match(Metadata_Mama_Fibro$ID, colnames(counts)  )), ]
  rownames(Metadata_Mama_Fibro) <- Metadata_Mama_Fibro$ID

  all(rownames(Metadata_Mama_Fibro) %in% colnames(counts))

  DGE <- DGEList(counts,
                 lib.size = colSums(counts),
                 norm.factors = calcNormFactors(counts),
                 samples = Metadata_Mama_Fibro,
                 genes = environ_genes)

  rownames(counts) <- environ_genes[,2]

  #write.xlsx(counts, file = sprintf("%s/Counts.xlsx", dirname(path_fc)))
  #write_rds(counts, file = sprintf("%s/Counts.rds", dirname(path_fc)))
  write.xlsx(DGE, file = sprintf("%s/DGEList.rds", dirname(path_fc)))
  return(DGE)
}


