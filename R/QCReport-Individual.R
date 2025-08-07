#' @title Generate QC report
#' @description kvjdfnkvjdf
#' @param patients_dir fskbfjsf
#' @export
#' @import rmarkdown
#' @import ggplot2
generateQCreport_Individual <- function(patient_dir) {

  file_list <- list.files(patient_dir, recursive = FALSE, full.names= TRUE)
  dirR1 <- file_list[endsWith(file_list, "R1_fastqc")]
  dirR2 <- file_list[endsWith(file_list, "R2_fastqc")]
  id <- basename(patient_dir)
  #if(file.exists(sprintf("%s/QCReport_Biota.pdf", patient_dir))) {
  #  return(message("The QC report of this patient has already been done!"))
  #}

  # Contenido del archivo Rmd
  rmd_content <- sprintf('
---
title: "Reporte de Calidad de Muestra"
output:
  pdf_document:
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(readr)
library(stringr)
create_FQCdata <- function(dat) {
  DBname <- read_lines(dat)
  len <- length(DBname)
  DBname <- DBname[2:len]
  DBname <- read.table(text = DBname, sep = "\t")
  colnames(DBname) <- strsplit(read_lines(dat)[2],"\t", fixed=TRUE)[[1]]
  colnames(DBname)[1] <- str_replace(colnames(DBname)[1], "#", "")
  return(DBname)
}

dir_fastqc_R1 <- "%s"
dir_fastqc_R2 <- "%s"
report_R1 <- read_file(paste0(dir_fastqc_R1, "/fastqc_data.txt", sep=""))
report <- as.data.frame(report_R1)
module_R1 <- str_split(report_R1, ">>")
#print(basename(dir_fastqc_R2))
report_R2 <- read_file(paste0(dir_fastqc_R2, "/fastqc_data.txt", sep=""))
module_R2 <- str_split(report_R2, ">>")


basic_stats <- create_FQCdata(module_R1[[1]][2])
Overrepresented_sequences <- create_FQCdata(module_R1[[1]][20])
cant_over <- nrow(Overrepresented_sequences)
basic_stats <- rbind(
  basic_stats,
  data.frame(
    Measure = "Overrepresented sequences",
    Value = cant_over,
    stringsAsFactors = FALSE
  )
)

#Per_tile_sequence_quality <- create_FQCdata(module_R1[[1]][6])
Per_sequence_quality_scores <- create_FQCdata(module_R1[[1]][8])
library(ggplot2)
plotPSQS <- ggplot(Per_sequence_quality_scores, aes(x = Quality, y = Count)) +
    geom_line() +
    labs(title = "Per Sequence Quality Scores",
         x = "Mean sequence quality",
         y = "") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_line(color = "white", size = 0.2),  # Líneas más claras y delgadas
        panel.grid.minor = element_blank(),  # Elimina las líneas menores
        panel.background = element_rect(fill = "white", color = NA)  # Fondo blanco sin bordes
      )

print(plotPSQS)



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

Per_base_sequence_quality$Median <- ((Per_base_sequence_quality_R1$Median + Per_base_sequence_quality_R2$Median) / 2)
Per_base_sequence_quality$`10th Percentile` <- ((Per_base_sequence_quality_R1$`10th Percentile`  + Per_base_sequence_quality_R2$`10th Percentile` ) / 2)
Per_base_sequence_quality$`90th Percentile` <- ((Per_base_sequence_quality_R1$`90th Percentile`  + Per_base_sequence_quality_R2$`90th Percentile` ) / 2)
Per_base_sequence_quality$`Lower Quartile` <- ((Per_base_sequence_quality_R1$`Lower Quartile`  + Per_base_sequence_quality_R2$`Lower Quartile` ) / 2)
Per_base_sequence_quality$`Upper Quartile` <- ((Per_base_sequence_quality_R1$`Upper Quartile`  + Per_base_sequence_quality_R2$`Upper Quartile` ) / 2)

plotPBSQ <- ggplot(Per_base_sequence_quality, aes(x = Base, y = Mean)) +
      geom_boxplot(width = 0.5, fill = "black", color = "black") +
      geom_point(aes(y = Median), color = "yellow", size = 2, position = position_dodge(width = 0.75)) +
      geom_errorbar(
        aes(ymin = `Lower Quartile`, ymax = `Upper Quartile`),
        width = 0.2,
        position = position_dodge(width = 0.75),
        color = "black"
      ) +
      geom_linerange(
        aes(ymin = `10th Percentile`, ymax = `90th Percentile`),
        position = position_dodge(width = 0.75),
        color = "black"
      ) +
      geom_line(aes(y = Mean, group = 1), position = position_dodge(width = 0.75), color = "blue", size = 1) +
      geom_hline(yintercept = 28, linetype = "dashed", color = "green", size = 0.7) +  # Add horizontal lines
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 0.7) +  # Add horizontal lines
      labs(title = "Per Base Sequence Quality",
           x = "Position in read (Base)",
           y = "Values") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_line(color = "white", size = 0.2),  # Líneas más claras y delgadas
        panel.grid.minor = element_blank(),  # Elimina las líneas menores
        panel.background = element_rect(fill = "white", color = NA)  # Fondo blanco sin bordes
      ) +  # Rotate x-axis labels
      scale_y_continuous(breaks = seq(0, 40, by = 2))+
      coord_cartesian(ylim = c(0, 41))  # Set y-axis limits

print(plotPBSQ)


```

```{r, echo=FALSE}
knitr::kable(basic_stats, caption = "Basic Statistics")
print(plotPBSQ)
print(plotPSQS)

```
  ', dirR1, dirR2)

rmd_file <- sprintf("%s/QCReport_Biota.Rmd", patient_dir)
writeLines(rmd_content, con = rmd_file)
#assign("plotPBSQ", plotPBSQ, envir = .GlobalEnv)
#assign("plotPSQS", plotPSQS, envir = .GlobalEnv)
#assign("basic_stats", basic_stats, envir = .GlobalEnv)

#sudo apt-get install texlive-full
#tinytex::install_tinytex()

id <- basename(patient_dir)
output_file <- sprintf("%s/QCReport_ID%s.pdf", patient_dir, id)
rmarkdown::render(rmd_file, output_file = output_file)
#rmarkdown::render(rmd_file)

}


plotsFastQC <- function(patients_dir, trimmed = FALSE, R= "R1R2") {

  #Separación de los modulos del fastqc para generar graficos
  dir_list <- list.dirs(path = patients_dir, full.names = FALSE, recursive = FALSE)
  cant_patients <- length(dir_list)

  list_modulo_PBSQ <- list() #Va a almacenar los dataframes para el grafico Per Base Sequence Quality

  #Recorro cada paciente para guardar las tablas para el grafico y sacar el promedio de las medias.
  for (p in 1:cant_patients) {

    #Toma los fastqc y se asegura de que estén unzippeados
    if (trimmed == FALSE) {

      file_list <- list.files(sprintf("%s/%s", patients_dir, dir_list[p]))

      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_001_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s", patients_dir, dir_list[p]))
      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_001_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s", patients_dir, dir_list[p]))

      file_list <- list.files(sprintf("%s/%s", patients_dir, dir_list[p]))

      if (!(length(nchar(file_list[endsWith(file_list, "R1_fastqc")])) == 0 )) {
        dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_fastqc")], sep="")
      } else if (!(length(nchar(file_list[endsWith(file_list, "R1_001_fastqc")])) == 0 )) {
        dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_001_fastqc")], sep="")
      }

      if(!(length(nchar(file_list[endsWith(file_list, "R2_fastqc")])) == 0 )) {
        dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_fastqc")], sep="")
      } else if(!(length(nchar(file_list[endsWith(file_list, "R2_001_fastqc")])) == 0 )) {
        dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_001_fastqc")], sep="")
      }

      dir_fastq_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_001.fastq.gz")], sep="")
      dir_fastq_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_001.fastq.gz")], sep="")


    } else {
      file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))

      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))
      file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fastqc.zip")], sep="")
      unzip(file_fastqc_zip, exdir = sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))

      file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))
      dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fastqc")], sep="")
      dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fastqc")], sep="")

      dir_fastq_R1 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fq.gz")], sep="")
      dir_fastq_R2 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fq.gz")], sep="")
    }

    if( length(nchar(dir_fastqc_R1)) == 0 | length(nchar(dir_fastqc_R2)) == 0 ) {
      stop("There is no fastQC file in this folder. Try trimmed == FALSE.")
    }
    print(basename(dir_fastqc_R1))


    library(readr)
    library(stringr)
    report_R1 <- read_file(paste0(dir_fastqc_R1, "/fastqc_data.txt", sep=""))
    report <- as.data.frame(report_R1)
    module_R1 <- str_split(report_R1, ">>")

    print(basename(dir_fastqc_R2))
    report_R2 <- read_file(paste0(dir_fastqc_R2, "/fastqc_data.txt", sep=""))
    module_R2 <- str_split(report_R2, ">>")

    basic_stats <- create_FQCdata(module_R1[[1]][2])
    #Valores de calidad:
    library(qckitfastq)
    #QC <- qual_score_per_read(dir_fastq_R1)
    #mean <- mean(QC$mu_per_read)
    #Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    #scores_qc[i, "Sample"] <- basename(p)
    #scores_qc[i, "Q_Mean_R1"] <- mean
    #scores_qc[i, "%_>=Q30_R1"] <- Q30

    #QC <- qual_score_per_read(dir_fastq_R2)
    #mean <- mean(QC$mu_per_read)
    #Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    #Per_tile_sequence_quality <- create_FQCdata(module_R1[[1]][6])
    Per_sequence_quality_scores <- create_FQCdata(module_R1[[1]][8])
    library(ggplot2)
    plotPSQS <- ggplot(Per_sequence_quality_scores, aes(x = Quality, y = Count)) +
      geom_line() +
      labs(title = "Per Sequence Quality Scores",
           x = "Mean sequence quality",
           y = "") +
      theme_minimal() +

      print(plotPSQS)

    #Per_base_sequence_content <- create_FQCdata(module_R1[[1]][10])
    #Per_sequence_GC_content <- create_FQCdata(module_R1[[1]][12])
    #Per_base_N_content <- create_FQCdata(module_R1[[1]][14])
    #Sequence_Length_Distribution <- create_FQCdata(module_R1[[1]][16])
    #Sequence_Duplication_Levels <- create_FQCdata(module_R1[[1]][18])
    Overrepresented_sequences <- create_FQCdata(module_R1[[1]][20])
    #Adapter_Content <- create_FQCdata(module_R1[[1]][22])

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

    Per_base_sequence_quality$Median <- ((Per_base_sequence_quality_R1$Median + Per_base_sequence_quality_R2$Median) / 2)
    Per_base_sequence_quality$`10th Percentile` <- ((Per_base_sequence_quality_R1$`10th Percentile`  + Per_base_sequence_quality_R2$`10th Percentile` ) / 2)
    Per_base_sequence_quality$`90th Percentile` <- ((Per_base_sequence_quality_R1$`90th Percentile`  + Per_base_sequence_quality_R2$`90th Percentile` ) / 2)
    Per_base_sequence_quality$`Lower Quartile` <- ((Per_base_sequence_quality_R1$`Lower Quartile`  + Per_base_sequence_quality_R2$`Lower Quartile` ) / 2)
    Per_base_sequence_quality$`Upper Quartile` <- ((Per_base_sequence_quality_R1$`Upper Quartile`  + Per_base_sequence_quality_R2$`Upper Quartile` ) / 2)

    plotPBSQ <- ggplot(Per_base_sequence_quality, aes(x = Base, y = Mean)) +
      geom_boxplot(width = 0.5, fill = "black", color = "black") +
      geom_point(aes(y = Median), color = "yellow", size = 1, position = position_dodge(width = 0.75)) +
      geom_errorbar(
        aes(ymin = `Lower Quartile`, ymax = `Upper Quartile`),
        width = 0.2,
        position = position_dodge(width = 0.75),
        color = "black"
      ) +
      geom_linerange(
        aes(ymin = `10th Percentile`, ymax = `90th Percentile`),
        position = position_dodge(width = 0.75),
        color = "black"
      ) +
      geom_line(aes(y = Mean, group = 1), position = position_dodge(width = 0.75), color = "black", size = 1) +
      geom_hline(yintercept = 28, linetype = "dashed", color = "green", size = 0.5) +  # Add horizontal lines
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", size = 0.5) +  # Add horizontal lines
      labs(title = "Per Base Sequence Quality",
           x = "Position in read (Base)",
           y = "Values") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +  # Rotate x-axis labels
      scale_y_continuous(breaks = seq(0, 40, by = 2))+
      coord_cartesian(ylim = c(0, 41))  # Set y-axis limits

    print(plotPBSQ)

    png(filename = sprintf("%s/%s/plotFastQC_PBSQ.png",  patients_dir, dir_list[p]), width = 800, height = 600)  # Ajusta el tamaño según tus necesidades
    print(plotPBSQ)
    dev.off()
    message("The plot from FastQC was saved")

    png(filename = sprintf("%s/%s/plotFastQC_PSQS.png", patients_dir, dir_list[p]), width = 800, height = 600)  # Ajusta el tamaño según tus necesidades
    print(plotPSQS)
    dev.off()
    message("The plot from FastQC was saved")

  }
}


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

