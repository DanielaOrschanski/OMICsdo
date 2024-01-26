
#Control de calidad para reporte: ########################################################3
runFastQC(patient_dir = "/media/16TBDisk/Daniela/MuestrasColon4/Co1525")
plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasColonGuada", trimmed = TRUE)
plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasColonGuada", trimmed = FALSE)

plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasMama", trimmed = TRUE)
plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasMama", trimmed = FALSE)


plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasColon", trimmed = TRUE)
plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasColon", trimmed = FALSE)

plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasColon2", trimmed = TRUE)
plotFastQC_PBSQ(patients_dir = "/media/16TBDisk/Daniela/MuestrasColon2", trimmed = FALSE)

remove.packages("TestLibrary")

counts_Colon <- FeatureCount_Report$counts



BiocManager::install("qckitfastq")
library(qckitfastq)

file_path <- "/media/16TBDisk/Daniela/MuestrasColon/Co1272/Co1272_combined_R1.fastq.gz"
q_per_read <- qual_score_per_read(file_path)

mu_per_read <- q_per_read$mu_per_read
percentage_q30_read <- mean(mu_per_read >= 30) * 100

mu_per_position <- q_per_read$mu_per_position
percentage_q30_position <- mean(mu_per_position >= 30) * 100

#################################################################################################
#Unir counts del FC:
count_mama_guada <- readRDS("/media/16TBDisk/Daniela/MuestrasMamaGuada/FeatureCount_Report.rds")
count_mama_guada <- count_mama_guada$counts
count_mama <- readRDS("/media/16TBDisk/Daniela/MuestrasMama/FeatureCount_Report.rds")
count_mama <- count_mama$counts
count_colon <- readRDS("/media/16TBDisk/Daniela/MuestrasColon/FeatureCount_Report_12Muestras.rds")
count_colon <- count_colon$counts
count_colon2 <- readRDS("/media/16TBDisk/Daniela/MuestrasColon2/FeatureCount_Report.rds")
count_colon2 <- count_colon2$counts
count_colon3 <- readRDS("/media/16TBDisk/Daniela/MuestrasColon3/FeatureCount_Report.rds")
count_colon3 <- count_colon3$counts
count_colon4 <- readRDS("/media/16TBDisk/Daniela/MuestrasColon4/FeatureCount_Report.rds")
count_colon4 <- count_colon4$counts
count_colon_guada <- readRDS("/media/16TBDisk/Daniela/MuestrasColonGuada/FeatureCount_Report.rds")
count_colon_guada <- count_colon_guada$counts

Todos_counts <- cbind(count_mama_guada, count_colon_guada, count_mama,
                      count_colon, count_colon2, count_colon3, count_colon4)

saveRDS(Todos_counts, file= "/media/16TBDisk/Daniela/Counts_91.rds")


#################################
# Comparacion de calidad
################################
BiocManager::install("qckitfastq")
library(qckitfastq)
library(readr)
library(stringr)

#Run FastQC
runFastQC(patient_dir = "/media/16TBDisk/Daniela/MuestrasColonGuada")

infile <- "/media/16TBDisk/Daniela/MuestrasColonGuada/Co0173/Co0173_combined_R1.fastq.gz"
QC <- qual_score_per_read(infile)
qual_score_per_read(infile)$q50_per_position[1:10]


#Calcular el mean manualmente:
patient_dir <- "~/Biota/MuestrasLeo/80"
PBSQ_Colombia <- ModuloPBSQ(patient_dir)
mean(PBSQ_Colombia$Mean)

#Recorro cada paciente para guardar las tablas para el grafico y sacar el promedio de las medias.
ModuloPBSQ <- function(patient_dir, trimmed= FALSE) {

  if (trimmed == FALSE) {
    file_list <- list.files(patient_dir)
    dir_fastqc_R1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "R1_001_fastqc")], sep="")
    dir_fastqc_R2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "R2_001_fastqc")], sep="")

  } else {
    file_list <- list.files(patient_dir)
    dir_fastqc_R1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1_fastqc")], sep="")
    dir_fastqc_R2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_2_fastqc")], sep="")
  }

  if( length(nchar(dir_fastqc_R1)) == 0 | length(nchar(dir_fastqc_R2)) == 0 ) {
    stop("There is no fastQC file in this folder. Try trimmed == FALSE.")
  }
  print(basename(dir_fastqc_R1))

  report_R1 <- read_file(paste0(dir_fastqc_R1, "/fastqc_data.txt", sep=""))
  module_R1 <- str_split(report_R1, ">>")


  print(basename(dir_fastqc_R2))
  report_R2 <- read_file(paste0(dir_fastqc_R2, "/fastqc_data.txt", sep=""))
  module_R2 <- str_split(report_R2, ">>")

  #Plot principal:
  Per_base_sequence_quality_R1 <- create_FQCdata(module_R1[[1]][4])
  order_fact <- Per_base_sequence_quality_R1$Base
  Per_base_sequence_quality_R1$Base <- factor(Per_base_sequence_quality_R1$Base, levels = order_fact)

  Per_base_sequence_quality_R2 <- create_FQCdata(module_R2[[1]][4])
  order_fact <- Per_base_sequence_quality_R2$Base
  Per_base_sequence_quality_R2$Base <- factor(Per_base_sequence_quality_R2$Base, levels = order_fact)

  #promedio la mean de R1 y R2 para poder tener una sola línea por cada paciente
  Per_base_sequence_quality <- Per_base_sequence_quality_R1
  Per_base_sequence_quality$Mean <- ((Per_base_sequence_quality_R1$Mean + Per_base_sequence_quality_R2$Mean) / 2)

  message(sprintf("The module of patient %s for PBSQ graph has been loaded", patient_dir))
  return(Per_base_sequence_quality)
}

mean(Per_base_sequence_quality$Mean)



#' @title create FastQC data
#' @description Prepare the data for plotting
#' @param dat data from module of FastQC output
create_FQCdata <- function(dat) {
  DBname <- read_lines(dat)
  len <- length(DBname)
  DBname <- DBname[2:len]
  DBname <- read.table(text = DBname, sep = "\t")
  colnames(DBname) <- strsplit(read_lines(dat)[2],"\t", fixed=TRUE)[[1]]
  colnames(DBname)[1] <- str_replace(colnames(DBname)[1], "#", "")
  return(DBname)
}
