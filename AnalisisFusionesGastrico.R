setwd("/media/4tb1/")
ena_urls <- read.delim("/media/4tb1/Daniela/Fusiones/Fusiones/Gastrico-80p/filereport_read_run_PRJNA508414_tsv.txt", header = FALSE)

colnames(ena_urls) <- ena_urls[1,]
ena_urls <- ena_urls[-1,]

ena_urls$fastq_split <- strsplit(ena_urls$fastq_ftp, ";")
ena_urls$fastq_1 <- sapply(ena_urls$fastq_split, `[`, 1)
ena_urls$fastq_2 <- sapply(ena_urls$fastq_split, `[`, 2)
fastq_list <- unlist(ena_urls[, c("fastq_1", "fastq_2")])
fastq_list <- unname(fastq_list)

path_dir <- "/media/4tb1/Daniela/Fusiones/Fusiones/Gastrico-80p"


dir.create(path_dir)
#fastq_list <- fastq_list[-c(10)]
i=162
for (i in 1:length(fastq_list)) {
  options(timeout = 1000000)
  fastq <- fastq_list[i]
  
  id <- strsplit(basename(fastq), split="_")[[1]][1]
  R <- ifelse(grepl(1, strsplit(basename(fastq), split="_")[[1]][2]), "R1", "R2")
  dir.create(sprintf("%s/%s", path_dir, id))
  
  file <- sprintf("%s/%s/%s_%s.fastq.gz", path_dir, id, id, R)
  if(file.exists(file)) {
    if(file.size(file) == 0) {
      tryCatch({
        download.file(fastq, destfile = sprintf("%s/%s/%s_%s.fastq.gz", path_dir, id, id, R))
        #cmd <- sprintf("wget -q '%s' -O %s/%s/%s_%s.fastq.gz", fastq, path_dir, id, id, R)
        #system(cmd)
      }, error = function(e) {
        message(sprintf("Error downloading file for %s_%s: %s", id, R, e$message))
      })
    }
  }
  
  
  if(!(file.exists(sprintf("%s/%s/%s_%s.fastq.gz", path_dir, id, id, R))) ) {
    
    tryCatch({
      download.file(fastq, destfile = sprintf("%s/%s/%s_%s.fastq.gz", path_dir, id, id, R))
      #cmd <- sprintf("wget -q '%s' -O %s/%s/%s_%s.fastq.gz", fastq, path_dir, id, id, R)
      #system(cmd)
    }, error = function(e) {
      message(sprintf("Error downloading file for %s_%s: %s", id, R, e$message))
    })
    
    #setwd(path_dir)
    
  } else {
    print(sprintf("La muestra %s %s ya fue descargada", id, R))
  }
  i <- i+1
}

path_dir = "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446"
path_dir = "/media/16TBDisk/Daniela/Fusiones/Gastrico-610"

RNAseqP(path_dir, fastQC_after_trim = FALSE, plot_FastQC_trim = FALSE, RunARRIBA = TRUE, RunFeatureCounts =  TRUE)


