#' @title index Reference genome with STAR
#' @description Index genome reference version HG38 with STAR.
#' @return Paths of FASTA and GTF(annotation) from the genome reference.
indexRefSTAR <- function(AnnotationHG38, FastaHG38) {
  message("ESTOY EN EL INDEX")
  soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  #soft_directory <- sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))

  #softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER') ))
  softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
  linea_software <- grep("(?i)HG38Index", softwares, ignore.case = TRUE, value = TRUE)

  if(length(nchar(linea_software)) == 0) {
    #Index genome reference
    message("The genome reference will be index with STAR, please be patient this process may take some minutes.")
    nThreads <- max(1,parallel::detectCores()-2)

    STAR <- downloadSTAR()

    # Cambiar permisos para lectura
    system(paste("chmod +r", FastaHG38), intern = TRUE)

    log <- system2(command = STAR,
                   args = c(paste0("--runThreadN ", nThreads),
                            "--runMode genomeGenerate",
                            paste0("--genomeDir " , sprintf("%s/HG38/index", soft_directory)),
                            paste0("--genomeFastaFiles ", FastaHG38),
                            "--sjdbOverhang 100 ",
                            "--genomeSAindexNbases 12",
                            paste0("--sjdbGTFfile ", AnnotationHG38)
                   ), stdout = TRUE)#out.file)

    #Writes down the paths in the txt
    HG38Index <- sprintf("%s/HG38/index", soft_directory)
    #HG38Index <- sprintf("%s/HG38/indexP", soft_directory)

    #softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
    softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
    softwares_actualizado <- c(softwares, sprintf("HG38Index %s", HG38Index))
    #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
    write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    return(HG38Index)

  } else {
    message("The reference genome was already indexed with STAR")
    HG38IndexSTAR <- sprintf("%s/HG38/index", soft_directory)
    return(HG38IndexSTAR)
  }
}


#' @title downloadHG38
#' @description Downloads the FASTA and the annotation of the genome reference version HG38.
#' @return Paths of FASTA and GTF(annotation) from the genome reference.
#' @import GEOquery
#' @export
downloadHG38 <- function() {
  message("ESTOY EN DOWNLOAD HG38")
  #soft_directory <- sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER'))
  soft_directory <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))


  #Checks if Annotation is already downloaded --------------------
  #softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
  softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

  linea_software <- grep("HG38Annotation", softwares, ignore.case = TRUE, value = TRUE)

  if(length(nchar(linea_software)) == 0) {
    #download Annotation
    message("HG38 annotation will be downloaded")
    dir.create(sprintf("%s/HG38", soft_directory))
    URL <- "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

    dir <- sprintf("%s/HG38", soft_directory)
    setwd(dir)
    system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
    AnnotationHG38 <- sprintf("%s/HG38/Homo_sapiens.GRCh38.110.gtf.gz", soft_directory)

    gunzip(AnnotationHG38, destname = gsub("[.]gz$", "", AnnotationHG38), overwrite = FALSE, remove = TRUE)
    AnnotationHG38 <<- sprintf("%s/HG38/Homo_sapiens.GRCh38.110.gtf", soft_directory)

    #Writes down the paths in the txt
    #softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
    softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    softwares_actualizado <- c(softwares, sprintf("HG38Annotation %s", AnnotationHG38))
    #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
    write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

  } else {
    message("The annotation for HG38 was already downloaded")
    AnnotationHG38 <<- strsplit(linea_software, " ")[[1]][[2]]
  }

  #Checks if FASTA is already downloaded --------------------------------------------------------
  linea_software <- grep("HG38FASTA", softwares, ignore.case = TRUE, value = TRUE)

  if(length(nchar(linea_software)) == 0) {
    message("HG38 FASTA will be downloaded")
    #download Fasta HG38
    #setwd(sprintf("%s/OMICsdoSof", Sys.getenv('R_LIBS_USER')))
    #URL <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz"

    URL <- "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
    system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)

    FastaHG38 <- sprintf("%s/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz", soft_directory)
    gunzip(FastaHG38, destname = gsub("[.]gz$", "", FastaHG38), overwrite = FALSE, remove = TRUE)
    FastaHG38 <<- sprintf("%s/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa", soft_directory)
    #FastaHG38 <<- sprintf("%s/HG38/GCA_000001405.15_GRCh38_genomic.fna", soft_directory)

    #Indexar con bwa la referencia:
    message("The reference will be indexed by BWA")
    system(sprintf("bwa index %s"), FastaHG38)
     #Writes down the paths in the txt
    #softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
    softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

    softwares_actualizado <- c(softwares, sprintf("HG38FASTA %s", FastaHG38))
    #softwares_actualizado <- c(softwares, sprintf("HG38FASTAP %s", FastaHG38))
    #write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", Sys.getenv('R_LIBS_USER')))
    write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

  } else {
    message("The fasta file for HG38 was already downloaded")
    FastaHG38 <<- strsplit(linea_software, " ")[[1]][[2]]
  }


  index_dir_STAR <- indexRefSTAR(AnnotationHG38, FastaHG38)

  return(c(AnnotationHG38, FastaHG38, index_dir_STAR))

}


#' @title downloadHG19
#' @description Downloads the FASTA and the annotation of the genome reference version HG38.
#' @return Paths of FASTA and GTF(annotation) from the genome reference.
#' @import GEOquery
#' @export
downloadHG19 <- function() {
  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))

  softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))
  linea_software <- grep("(?i)HG19FASTA", softwares, ignore.case = TRUE, value = TRUE)

  if(length(nchar(linea_software)) == 0) {

    message("The genome reference hg19 will be downloaded and indexed, please be patient this process may take some minutes.")

    dir.create(sprintf("%s/HG19", omicsdo_sof))
    URL <- "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
    dir <- sprintf("%s/HG19", omicsdo_sof)
    system2("wget", args = c(URL, "-P", dir), wait = TRUE, stdout = NULL, stderr = NULL)
    system2("gzip" , sprintf("-d %s/%s", dir, basename(URL)))

    reference <- substr(basename(URL), 1, nchar(basename(URL)) - 3)
    referenceN <- stringr::str_replace(reference, "fa", "fasta")
    referenceN <- sprintf("%s/%s", dir, referenceN)
    reference <- sprintf("%s/%s", dir, reference)

    system2("mv", sprintf("%s %s", reference, referenceN))

    BWA <- downloadBWA()
    GATK <- downloadGATK()
    SAMTOOLS <- downloadSamtools()

    #system2(BWA, sprintf("index -a is %s", referenceN))
    message("The reference will be indexed by BWA")
    system(sprintf("bwa index %s"), referenceN)

    #indexarla con STAR tambien

    system2("java", sprintf("-jar %s CreateSequenceDictionary -R %s", GATK, referenceN))
    system2(SAMTOOLS, sprintf("faidx %s", referenceN))

    FastaHG19ADN <<- referenceN
    softwares_actualizado <- c(softwares, sprintf("HG19FASTA %s", FastaHG19ADN))
    write(softwares_actualizado, file = sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

  } else {
    message("The reference genome hg19 was already indexed with BWA")
    FastaHG19ADN <<- linea_software
  }
  data("exons.hg19")
  bed <<- exons.hg19
  return(FastaHG19ADN)

}



