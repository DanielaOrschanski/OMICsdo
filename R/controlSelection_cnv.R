#' @import ExomeDepth
#' @title Control Selection for CNVs
#' @description Using the ExomeDepth algorithm to select the patient that is more likely to be "normal" (no CNVs).
#' @param id number or character that indicates which patient will be analyzed.
#' It has to be the same as its id in the data base.
#' @param PatientsDB database in dataframe format that contains the counts of
#' the genes of any patient (normal or not) that has been analyzed.
#' @param bed is a dataframe that has start, end and name of exons within the panel selected.
#' @return controles dataframe that indicates which samples has been selected as references and it frequency.
#' @export

controlSelection_cnv <- function(PatientsDB, bed) {

  controles <- data.frame(matrix(NA, nrow = length(colnames(PatientsDB)[-c(1:4)]), ncol = 2))

  for (i in 1:length(colnames(PatientsDB)[-c(1:4)])) {
    #id <- "BM22-46334"
    id <- colnames(PatientsDB)[-c(1:4)][i]
    print(id)
    indice <- which(colnames(PatientsDB) == id)
    my.test <- as.numeric(PatientsDB[, indice, drop = TRUE])
    my.ref.samples <- colnames(PatientsDB)[-c(1,2,3,4, indice)]
    my.reference.selected <- apply(X = PatientsDB[, my.ref.samples, drop = FALSE],
                                   MAR = 1,
                                   FUN = sum)

    my.reference.set <- as.matrix(PatientsDB[, my.ref.samples])
    my.choice <- select.reference.set (test.counts = my.test,
                                       reference.counts = my.reference.set,
                                       bin.length = (PatientsDB$end - PatientsDB$start)/1000,
                                       n.bins.reduced = 10000)

    print(my.choice[[1]])
    controles[i, 1] <- id
    controles[i, 2] <- paste(my.choice[[1]], collapse=",")

  }
  colnames(controles) <- c("ID", "Referencias")
  referencias <- unlist(strsplit(controles$Referencias, ","))

  frecuencia_IDs <- as.data.frame(table(referencias))
  colnames(frecuencia_IDs)[1] <- "ID"
  controles <- merge(controles, frecuencia_IDs, by="ID")

  return(controles)
}
