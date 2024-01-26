#' @title downloadHG38
#' @description Downloads the FASTA and the annotation of the genome reference version HG38.
#' @return Paths of FASTA and GTF(annotation) from the genome reference.
#' @import GEOquery
#' @export

runMIXTURE <- function() {

  data(LM22)
  data("TIL10")
  LM22 <- LM22
  TIL10 <- TIL10

  rownames(LM22)[1:10]

  # Instalar y cargar la base de datos org.Hs.eg.db
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
  }
  library(org.Hs.eg.db)
  BiocManager::install("annotate")
  library(annotate)
  #https://rdrr.io/bioc/GeneAnswers/man/getSymbols.html


  # Obtener el símbolo del gen para los primeros 10 identificadores
  gene_ids <- c("100287102", "653635", "102466751", "100302278", "645520", "79501", "729737", "102725121", "102723897", "102465909")

  # Consultar la base de datos para obtener los símbolos de genes
  gene_symbols <- mapIds(org.Hs.eg.db, gene_ids, "SYMBOL", "ENTREZID")

  # Imprimir los resultados
  print(gene_symbols)



  mix.test <- MIXTURE(expressionMatrix = LM22,          #N x ncol(signatureMatrix) gene expresion matrix to evaluate
                      ##rownames(M) should be the GeneSymbols
                      signatureMatrix = LM22,                 #the gene signature matrix (W) such that M = W*betas'
                      #(i.e the LM22 from Newman et al)
                      iter = 1000,                            #iterations for the statistical test (null distribution)
                      functionMixture = nu.svm.robust.RFE,    #cibersort, nu.svm.robust.rfe, ls.rfe.abbas,
                      useCores = 10L,                         #cores for parallel processing/ if using windows set to 1
                      verbose = TRUE,                         #TRUE or FALSE messages
                      nullDist = "PopulationBased",           #"none" or "PopulationBased" if the statistical test should
                      #be performed
                      fileSave = "MIXTURE_FILE_LM22.xlsx")    #EXCEL file name to store the results

  save(mix.test, file = "MIXTURE_FILE_LM22.RData") #save full list as an RData object.

  head(GetMixture(exp.mix))
}



