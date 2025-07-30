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
setwd(path_dir)
#fastq_list <- fastq_list[-c(10)]
i=1
for (i in 1:length(fastq_list)) {
  options(timeout = 100000000000)
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
path_dir = "/media/4tb1/Daniela/Fusiones/Fusiones/Gastrico-80p"
length(list.dirs(path_dir))

library(OMICsdo)
RNAseqP(path_dir, fastQC_after_trim = FALSE, plot_FastQC_trim = FALSE, RunARRIBA = TRUE, RunFeatureCounts =  TRUE)
library(readr)
Metadata <- read_csv("/media/16TBDisk/Daniela/Fusiones/Gastrico-610/Metadata-Gastrico-610.csv")
#Metadata <- read_csv("/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SraRunTable.csv")

colnames(Metadata)[1] <- "ID"
Metadata_sin_celllines <- Metadata[-which(Metadata$source_name == "Gastric cancer cell line"),]

out <- fusionStats(
  #patients_dir = "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446",
  patients_dir = "/media/16TBDisk/Daniela/Fusiones/Gastrico-610",
  #Metadata = Metadata,
  Metadata = Metadata_sin_celllines,
  #Metadata = Metadata_apareados,
  #Apareados = TRUE,
  #group = "disease_stage",
  #group = "Subtipo",
  #group = "Histology",
  group = "tissuetype",
  #group = "Pathcode_min=0_wide=1",
  #cohorte = "Tiroides-446",
  cohorte = "Gastrico-610",
  sobrevida = FALSE)

table(Metadata_sin_celllines$tissuetype)

library(readxl)
Todos_FusionReports_Gastrico_610 <- read_excel("/media/16TBDisk/Daniela/Fusiones/Gastrico-610/Todos-FusionReports_Gastrico-610.xlsx")
colnames(Todos_FusionReports_Gastrico_610)[1] <- "ID"
Fusions_Metadata <- merge(Todos_FusionReports_Gastrico_610, Metadata, by = "ID")
Fusions_Metadata <- Fusions_Metadata[which(Fusions_Metadata$source_name == "Gastric cancer cell line"), ]
#CCHCR1::HLA-L
Fusions_Metadata <- Fusions_Metadata[which(Fusions_Metadata$gene1 == "CCHCR1" & Fusions_Metadata$gene2 == "HLA-L"),]
Fusions_Metadata <- Fusions_Metadata[-which(duplicated(Fusions_Metadata$ID)),]
table(Fusions_Metadata$ID)
table(Fusions_Metadata$Histology)
table(Fusions_Metadata$tissuetype)
table(Fusions_Metadata$source_name)

##########################################################
# TIROIDS
Metadata <- read_csv("/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SraRunTable.csv")
Metadata <- Metadata[Metadata$sample_type == "RNASeq",]

Metadata$disease_stage
colnames(Metadata)[colnames(Metadata) == "OS(0=_dead"] <- "tiempo"
Metadata$tiempo <- as.numeric(Metadata$tiempo)
#Trunco el tiempo a 50 para ver si me da significativo:
#Metadata$tiempo <- ifelse(Metadata$tiempo > 50, 50, Metadata$tiempo)

colnames(Metadata)[colnames(Metadata) == "_1_=_alive)" ] <- "evento"
Metadata$evento <- ifelse(Metadata$evento == "No", 0, ifelse(Metadata$evento == "Yes", 1, Metadata$evento))
Metadata$evento <- as.factor(Metadata$evento)


colnames(Metadata)[1] <- "ID"
Apareados <- Metadata
InfoImp <- Apareados$`Library Name`
Apareados$InfoImp <- InfoImp
Apareados$Subtipo <- sapply(strsplit(InfoImp, split = "_"), `[`, 1)
Apareados$NuevoID <- sapply(strsplit(InfoImp, split = "_"), `[`, 2)
Apareados$Tejido <- sapply(strsplit(InfoImp, split = "_"), `[`, 3)

colnames(Apareados)
Apareados <- Apareados[, c("ID","InfoImp" , "NuevoID", "Subtipo", "Tejido", "sex", "age", "disease", "disease_stage", "pNstage", "pT_stage" )]

Metadata$InfoImp <- Apareados$InfoImp
Metadata$Subtipo <- Apareados$Subtipo
Metadata$NuevoID <- Apareados$NuevoID
Metadata$Tejido <- Apareados$Tejido

#Los que tienen mismo Subtipo y mismo Nuevo ID son una misma muestra apareada!!!!:

Metadata[which(Metadata$`Pathcode_min=0_wide=1` == "missing"),] <- NA
Metadata_apareados <- Metadata[which(Metadata$NuevoID %in% c("1", "10", "15", "18", "25")),]
Metadata_apareados <- Metadata_apareados[order(Metadata_apareados$NuevoID), ]
Metadata_apareados <- Metadata_apareados[-c(2,3,5,8,11),]


Metadata_tumorales <- Metadata[which(Metadata$Tejido == "RNAtumor"), ]
Metadata_tumorales$tissue
Metadata$tissue

Metadata$tiempo
out <- fusionStats(
  patients_dir = "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446",
  #Metadata = Metadata,
  Metadata = Metadata_tumorales,
  #Metadata = Metadata_apareados,
  #Apareados = TRUE,
  group = "Subtipo",
  cohorte = "Tiroides-446",
  sobrevida = TRUE)
