#' @title download arcasHLA
#' @description Executes the
#' @export
download_arcasHLA <- function() {
  install_github("https://github.com/RabadanLab/arcasHLA")

  devtools::install_github("RabadanLab/arcasHLA")
  remotes::install_github("RabadanLab/arcasHLA")

  #samtools sort $DIR/$item/$item.bam -o $DIR/$item/${item}_sorted.bam
  #arcasHLA extract $DIR/$item/${item}_sorted.bam -o $DIR/$item
  #arcasHLA genotype $DIR/$item/${item}_sorted.extracted.1.fq.gz $DIR/$item/${item}_sorted.extracted.2.fq.gz -o $DIR/$item

  install.packages("reticulate")
  library(reticulate)

  Sys.which("python")
  use_python("/home/daniela/miniconda3/bin/python")
  py_install("git+https://github.com/RabadanLab/arcasHLA.git")
  arcasHLA <- import("arcasHLA")
}

