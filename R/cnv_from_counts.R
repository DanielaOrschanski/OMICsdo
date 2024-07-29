#' @import ExomeDepth
#' @title Get CNV calls from gene counts
#' @description Using the ExomeDepth algorithm to detect deletions or duplications
#' of regions within the genes of the patient by comparing it against the normal
#' patients (my.ref.samples).
#' @param id number or character that indicates which patient will be analyzed.
#' It has to be the same as its id in the data base.
#' @param bed bedfile that has information about the genes of the mitochondrial
#' genome. It's first 3 columns must be: chromosome, start, end. This file is
#' given on this package as bedfileMito.RDS in the Data folder.
#' @param ControlDB database in dataframe format that contains the counts of the
#' genes of normal patients that are used as a reference for making the analysis.
#' @param PatientsDB database in dataframe format that contains the counts of
#' the genes of any patient (normal or not) that has been analyzed.
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain
#' from the normal copy number state to either a deletion or a duplication.
#' The default (0.0001) expect approximately 20 CNVs genome-wide.
#' @return CNV calls (dataframe that indicated the CNVs detected).
#' @examples
#' CNV_calls_102 <- CNVfromCounts(102, bed_file, ControlDB, PatientsDB, minoverlap=0.001, transit = 0,7 )
#' @export

cnv_from_counts <- function(id, bed, ControlDB, PatientsDB, minoverlap, transit) {

  indice <- which(colnames(PatientsDB) == id)
  my.test <- as.numeric(PatientsDB[, indice, drop = TRUE])

  my.ref.samples <- colnames(ControlDB)[-c(1:4)]
  my.reference.selected <- apply(X = ControlDB[, my.ref.samples, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)

  #Creates a ExomeDepth object with test and reference:
  all.exons <- new('ExomeDepth',
                   test = my.test ,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons, #ExomeDepth object
                        transition.probability = transit,
                        chromosome = bed$chromosome,
                        start = bed$start,
                        end = bed$end,
                        name = bed$name)

  exons.GRanges <- GenomicRanges::GRanges(seqnames = bed$chromosome,
                                          IRanges::IRanges(start = bed$start, end = bed$end),
                                          names = bed$name)

  #Add names of genes within each variation:
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.GRanges,
                             #reference.annotation = bed,
                             min.overlap = minoverlap, #mientras menor es, mas genes
                             column.name = 'name.gene')

  CNV_calls <- all.exons@CNV.calls

  #Order by size of variation:
  CNV_calls$size <- abs(CNV_calls$start - CNV_calls$end)
  CNV_calls <- CNV_calls[order(CNV_calls$size, decreasing = TRUE),]
  CNV_calls_original <- CNV_calls
  CNV_calls <- CNV_calls[-which(CNV_calls$BF < 0),]
  #CNV_calls_f <- CNV_calls[-which(is.na(CNV_calls$name.gene)),]

  #Agrego anotacion de genes ---------------------------
  CNV_calls$Exons <- NA
  CNV_calls$Genes <- NA
  i=1
  for (i in 1:nrow(CNV_calls)) {
    st <- CNV_calls$start[i]
    en <- CNV_calls$end[i]
    chr <- CNV_calls$chromosome[i]
    exones_afectados <- bed$name[which(bed$chromosome == chr & bed$start >= st & bed$end <= en)]
    CNV_calls$Exons[i] <- paste(exones_afectados, collapse = ", ")
    #length(exones_afectados) == CNV_calls$nexons[i]
    CNV_calls$nexons[i] <- length(exones_afectados)

    genes_afectados <- strsplit(exones_afectados, split="_")
    genes_afectados <- lapply(genes_afectados, function(x) x[1])
    genes_afectados <- unlist(unique(genes_afectados))
    CNV_calls$Genes[i] <- paste(genes_afectados, collapse = ", ")
    CNV_calls$ngenes[i] <- length(genes_afectados)
  }

  CNV_calls <- CNV_calls[, -which(colnames(CNV_calls)== "name.gene")]


  Annot <- all.exons@annotations


  return(list(CNV_calls, all.exons))

}
