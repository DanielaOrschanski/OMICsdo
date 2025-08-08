# REPORTE DE CALIDAD GRUPAL

#' @title generateCohortQCreport
#' @description kvjdfnkvjdf
#' @param patients_dir fskbfjsf
#' @export
#' @import rmarkdown
#' @import ggplot2

#generateCohortQCreport(patients_dir <- "/media/8tb02/Daniela/Leucemia",
#                       bbmap = TRUE,
#                       path_bbmap = "/media/8tb02/Daniela/Leucemia/bbmap",
#                       index_rsubread = "/media/8tb02/Daniela/Leucemia/Mycoplasma/indexRsubread/mycoplasma_index",
#                       rate = "0.6",
#                       rmd_content = rmd_content)


generateCohortQCreport <- function(patients_dir, bbmap, path_bbmap, index_rsubread, rate, rmd_content) {

  # Contenido del archivo Rmd
  rmd_content <- rmd_content
  rmd_file <- sprintf("%s/FinalQCReport.Rmd", patients_dir)
  writeLines(rmd_content, con = rmd_file)

  #sudo apt-get install texlive-full
  output_file <- sprintf("%s/FinalQCReport.pdf", patients_dir)

  rmarkdown::render(rmd_file, output_file = output_file)
  #rmarkdown::render(rmd_file)

}



runBBmap <- function(path_bbmap, patient_dir, rate) {
  id <- basename(patient_dir)
  print(id)
  file_list <- list.files(patient_dir, recursive= FALSE, full.names = TRUE)
  R1 <- file_list[endsWith(file_list, "R1.fastq.gz")]
  R2 <- file_list[endsWith(file_list, "R2.fastq.gz")]

  bbmap_dir <- sprintf("%s/AlineadoMycoplasmaRsubread", patient_dir)
  dir.create(bbmap_dir)
  R1_bbmap <- sprintf("%s/AlineadoMycoplasmaRsubread/%s_bbmap_%s_R1.fastq.gz", patient_dir, id, rate)
  R2_bbmap <- sprintf("%s/AlineadoMycoplasmaRsubread/%s_bbmap_%s_R2.fastq.gz", patient_dir, id, rate)

  if(!file.exists(R1_bbmap) | !file.exists(R2_bbmap)) {
    system(sprintf(
      "%s/reformat.sh in1=%s in2=%s out1=%s out2=%s samplerate=%s",
      path_bbmap, R1, R2, R1_bbmap, R2_bbmap, rate ))

    #system(sprintf("%s/reformat.sh in=%s out=%s samplerate=%s", path_bbmap, R2, R2_bbmap, rate))
  } else {
    message("BBmap files have already been generated")
  }

  #reformat.sh in= 5_S36_L001_R2_001.fastq.gz out= 5_S15_L001_R2_001.fastq.gz samplerate=0.6
  return(list(R1_bbmap, R2_bbmap))
}

mapeoMycoplasma_Grupal <- function(dir_patients, index_rsubread, bbmap = FALSE, path_bbmap = NA, rate = "0.6") {

  patients <- list.dirs(dir_patients, full.names = FALSE, recursive = FALSE)
  mapeos_mycoplasma <- data.frame(ID = rep(0, length(patients)), Mapeo = rep(0,length(patients)))
  #p=2
  for (p in 1:length(patients)) {
    id <- patients[p]
    print(id)
    patient_dir <- paste(dir_patients, "/", id, sep ="")

    out <- mapeoMycoplasma_Individual(patient_dir = patient_dir,
                                      path_bbmap = path_bbmap,
                                      index_rsubread = index_rsubread,
                                      bbmap = bbmap, rate = rate)

    mapeos_mycoplasma$Mapeo[p] <- out$Mapeo
    mapeos_mycoplasma$ID[p] <- out$ID

  }
  return(mapeos_mycoplasma)
}


mapeoMycoplasma_Individual <- function(patient_dir, index_rsubread, bbmap = FALSE, path_bbmap = NA, rate = "0.6") {

  id <- basename(patient_dir)
  print(id)

  if(bbmap == FALSE) { #Si no se habilita el bbmap, el mycofree se hace sobre los fastq originales --> Tarda mas tiempo!
    file_list <- list.files(patient_dir, recursive= FALSE, full.names = TRUE)
    R1 <- file_list[endsWith(file_list, "R1.fastq.gz")]
    R2 <- file_list[endsWith(file_list, "R2.fastq.gz")]
    message("Mycofree va a ejecutarse sobre los files originales")
    rate <- "100"

  } else {
    out <- runBBmap(path_bbmap = path_bbmap, patient_dir = patient_dir, rate = rate)
    R1 <- out[[1]]
    R2 <- out[[2]]
    message(sprintf("Mycofree va a ejecutarse sobre los files BBmap con rate %s", rate))
  }

  #R1 <- sprintf("%s/%s/%s_R1.fastq.gz", dir_patients, patients[p], patients[p])
  #R2 <- sprintf("%s/%s/%s_R2.fastq.gz", dir_patients, patients[p], patients[p])

  path <- sprintf("%s/AlineadoMycoplasmaRsubread", patient_dir)

  library(Rsubread)
  if (file.exists(paste(path, "/", id, "_rate", rate ,"_alignRsubread.BAM.summary", sep=""))) {
    print("ya se le hizo el alineamiento")
  } else {
    print("no se le hizo el alinemiento")
    dir.create(path)
    align(
      index = index_rsubread, #Tiene que tener el nombre de todos los archivos del index, no de la carpeta nomas.
      readfile1 = R1,
      readfile2 = R2,
      input_format = "FASTQ",
      output_file = paste(path, "/",id, "_rate", rate ,"_alignRsubread.BAM", sep=""),
      type = "rna",
      nthreads = 15
    )
  }
  summary <- read.delim(paste(path, "/", id, "_rate", rate,"_alignRsubread.BAM.summary", sep=""), header = FALSE)
  total <- summary[1,2]
  mapped <- summary[3,2]
  mapeo_unico <- round(mapped/total*100, 2)

  mapeo_mycoplasma <- data.frame(ID = id, Mapeo = mapeo_unico)

  return(mapeo_mycoplasma)
}

create_FQCdata <- function(dat) {
  DBname <- read_lines(dat)
  len <- length(DBname)
  DBname <- DBname[2:len]
  DBname <- read.table(text = DBname, sep = "\t")
  colnames(DBname) <- strsplit(read_lines(dat)[2],"\t", fixed=TRUE)[[1]]
  colnames(DBname)[1] <- str_replace(colnames(DBname)[1], "#", "")
  return(DBname)
}

calculateQCMetricsSTAR <- function(patients_dir, trimmed = FALSE) {

  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)
  metricasSTAR <- data.frame("Categoria" = c(), "Valor" = c())

  for (p in 1:cant_patients) {
    id <- basename(dir_list[p])

    if(trimmed == TRUE) {
      log_final <- sprintf("%s/trimmed/%sLog.final.out", dir_list[p], id)
    } else{
      log_final <- sprintf("%s/%sLog.final.out", dir_list[p], id)
    }

    log_file <- readLines(log_final)
    df <- data.frame()
    #Genero un df por cada paciente
    for (line in log_file) {
      parts <- strsplit(line, "\\s*\\|\\s*")[[1]]

      if (length(parts) == 2) {
        # Extraer el nombre de la columna (izquierda del "|") y el valor de la fila (derecha del "|")
        column_name <- trimws(parts[1])
        value <- trimws(parts[2])
        df <- rbind(df, c(column_name, value))
      }
    }

    #Acumulo esos df uno por cada columna
    colnames(df) <- c("Categoria", "Valor")

    #Si es la primera que se hace que se guarden ambas columnas, sino solo la que tiene valores
    ifelse(p==1, metricasSTAR <- df, metricasSTAR <- cbind(metricasSTAR, df$Valor))

    colnames(metricasSTAR)[p+1] <- id
  }

  metricasST <- as.data.frame(t(metricasSTAR[5:9, ]))
  colnames(metricasST) <- metricasST[1,]
  metricasST <- metricasST[-1,]
  metricasST <- cbind(Sample = rownames(metricasST), metricasST)
  colnames(metricasST)[1] <- "ID"
  write.xlsx(metricasST, file = sprintf("%s/MetricasQCSTAR.xlsx", patients_dir))

  metrics_long <- mutate(
    metricasST,
    across(
      .cols = -ID,
      .fns = ~ as.numeric(gsub("%", "", .))
    )
  )

  #return(metricasST)
  return(metrics_long)
}

plotFastQC_PBSQ <- function(patients_dir, trimmed = FALSE, R= "R1R2") {

  #Separación de los modulos del fastqc para generar graficos
  dir_list <- list.dirs(path = patients_dir, full.names = FALSE, recursive = FALSE)
  cant_patients <- length(dir_list)

  list_modulo_PBSQ <- list() #Va a almacenar los dataframes para el grafico Per Base Sequence Quality

  #Recorro cada paciente para guardar las tablas para el grafico y sacar el promedio de las medias.
  #p=1
  for (p in 1:cant_patients) {

    #Toma los fastqc y se asegura de que estén unzippeados
    if (trimmed == FALSE) {
      file_list <- list.files(sprintf("%s/%s", patients_dir, dir_list[p]))

      #file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1_fastqc.zip")], sep="")
      #unzip(file_fastqc_zip, exdir = sprintf("%s/%s", patients_dir, dir_list[p]))
      #file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2_fastqc.zip")], sep="")
      #unzip(file_fastqc_zip, exdir = sprintf("%s/%s", patients_dir, dir_list[p]))

      #file_list <- list.files(sprintf("%s/%s", patients_dir, dir_list[p]))

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

      dir_fastq_R1 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R1.fastq.gz")], sep="")
      dir_fastq_R2 <- paste0(patients_dir, "/", dir_list[p],"/", file_list[endsWith(file_list, "R2.fastq.gz")], sep="")


    } else {
      file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))

      #file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fastqc.zip")], sep="")
      #unzip(file_fastqc_zip, exdir = sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))
      #file_fastqc_zip <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fastqc.zip")], sep="")
      #unzip(file_fastqc_zip, exdir = sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))

      #file_list <- list.files(sprintf("%s/%s/trimmed", patients_dir, dir_list[p]))
      dir_fastqc_R1 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fastqc")], sep="")
      dir_fastqc_R2 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fastqc")], sep="")

      dir_fastq_R1 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_1_fq.gz")], sep="")
      dir_fastq_R2 <- paste0(patients_dir, "/", dir_list[p],"/trimmed/", file_list[endsWith(file_list, "val_2_fq.gz")], sep="")
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
    #Per_base_sequence_quality$Mean <- ((Per_base_sequence_quality_R1$Mean + Per_base_sequence_quality_R2$Mean) / 2)

    if (R == "R2") {
      Per_base_sequence_quality$Mean <- Per_base_sequence_quality_R2$Mean
    } else if ( R == "both" | R == "R1R2") {
      Per_base_sequence_quality$Mean <- ((Per_base_sequence_quality_R1$Mean + Per_base_sequence_quality_R2$Mean) / 2)
    } else if (R == "R1") {
      Per_base_sequence_quality$Mean <- Per_base_sequence_quality_R1$Mean
    }

    Per_base_sequence_quality$Median <- ((Per_base_sequence_quality_R1$Median + Per_base_sequence_quality_R2$Median) / 2)
    Per_base_sequence_quality$`10th Percentile` <- ((Per_base_sequence_quality_R1$`10th Percentile`  + Per_base_sequence_quality_R2$`10th Percentile` ) / 2)
    Per_base_sequence_quality$`90th Percentile` <- ((Per_base_sequence_quality_R1$`90th Percentile`  + Per_base_sequence_quality_R2$`90th Percentile` ) / 2)
    Per_base_sequence_quality$`Lower Quartile` <- ((Per_base_sequence_quality_R1$`Lower Quartile`  + Per_base_sequence_quality_R2$`Lower Quartile` ) / 2)
    Per_base_sequence_quality$`Upper Quartile` <- ((Per_base_sequence_quality_R1$`Upper Quartile`  + Per_base_sequence_quality_R2$`Upper Quartile` ) / 2)

    list_modulo_PBSQ[[p]] <- Per_base_sequence_quality
    message(sprintf("The module of patient number %s for PBSQ graph has been loaded", p))


    #Valores de calidad:
    #file_path <- sprintf("%s/%s", patients_dir, dir_list[p] )
    #library(qckitfastq)
    #q_per_read <- qual_score_per_read(file_path)

    #mu_per_read <- q_per_read$mu_per_read
    #percentage_q30_read <- mean(mu_per_read >= 30) * 100

    #mu_per_position <- q_per_read$mu_per_position
    #percentage_q30_position <- mean(mu_per_position >= 30) * 100

  }

  #Ploteo una linea por paciente superpuestas en un mismo grafico -----------------
  #Base del plot
  plot <- ggplot(Per_base_sequence_quality, aes(x = Base, y = Mean)) +
    geom_boxplot(width = 0.5, fill = "lightblue", color = "black", alpha= 0.5) +
    geom_point(aes(y = Median), color = "blue", size = 1, position = position_dodge(width = 0.75)) +
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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(),  # Quita las líneas principales de la cuadrícula
          panel.grid.minor = element_blank()   # Quita las líneas menores de la cuadrícula
    ) +  # Rotate x-axis labels
    scale_y_continuous(breaks = seq(0, 40, by = 2))+
    coord_cartesian(ylim = c(0, 41))  # Set y-axis limits

  #Genera un dataframe para poder agregar una linea por cada paciente:
  combined_dataframe <- do.call(cbind, lapply(seq_along(list_modulo_PBSQ), function(i) {
    df <- list_modulo_PBSQ[[i]]
    names(df) <- paste(names(df), "_", dir_list[i], sep = "_")
    df
  }))

  selected_columns <- combined_dataframe[, grepl("^Mean", names(combined_dataframe))]
  selected_columns$Base <-combined_dataframe[,1]
  df.long <- melt(selected_columns,id.vars="Base")

  # Crea el gráfico completo con una linea por paciente:
  #plot2 <- plot +
  #  geom_line(data = df.long, aes(x = Base, y = value, color = variable, group = variable), size = 1) +
  #  labs(title = "Per Base Sequence Quality", x = "Base", y = "Quality Value")

  plot2 <- plot +
    geom_line(data = df.long, aes(x = Base, y = value, color = variable, group = variable), size = 1) +
    labs(title = "Per Base Sequence Quality", x = "Base", y = "Quality Value") +
    guides(color = guide_legend(title = NULL)) +  # Quita el título de la leyenda
    theme(
      plot.title = element_text(size = 8),  # Tamaño de título (aunque está nulo)
      legend.text = element_text(size = 6),    # Tamaño del texto de las anotaciones
      axis.text.x = element_text(size = 7)
    )

  indicador <- ifelse(trimmed == TRUE, "_trimmeado", "_SINtrimmeado")
  #png(filename = sprintf("%s/plotFastQC_PBSQ%s.png", patients_dir, indicador), width = 800, height = 600)  # Ajusta el tamaño según tus necesidades
  #print(plot2)
  #dev.off()
  return(plot2)
}

calculateQCMetricsSeq <- function(patients_dir, trimmed = FALSE){
  library(qckitfastq)
  trimeado <- ifelse(trimmed == TRUE, "trimmed", "")
  file= paste(patients_dir, "/MetricasQC_", trimeado,".xlsx", sep ="")

  if(file.exists(file)) {
    metrics_seq <- read_excel(file)
    return(metrics_seq)
  }

  scores_qc <- data.frame("Sample"= c(), "Q_Mean_R1" = c(), "%_>=Q30_R1"= c(),
                          "Q_Mean_R2" = c(), "%_>=Q30_R2"= c())
  i=1
  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  cant_patients <- length(dir_list)

  for (p in dir_list) {
    #p <- dir_list[[i]]
    print(p)
    if  (trimmed == TRUE) {
      file_list <- list.files(sprintf("%s/trimmed/", p))
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "val_1.fq.gz")])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("val_1.fq%s", gzip))], sep="")
      fileR2 <- paste0(p, "/trimmed/", file_list[endsWith(file_list, sprintf("val_2.fq%s", gzip))], sep="")
    } else {
      file_list <- list.files(p)
      gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
      fileR1 <- paste0(p, "/", file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))], sep="")
      fileR2 <- paste0(p, "/", file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))], sep="")
    }

    if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
      stop("There are no fastq files in this directory")
    }

    QC <- qual_score_per_read(fileR1)
    mean <- mean(QC$mu_per_read)
    Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    scores_qc[i, "Sample"] <- basename(p)
    scores_qc[i, "Q_Mean_R1"] <- mean
    scores_qc[i, "%_>=Q30_R1"] <- Q30

    QC <- qual_score_per_read(fileR2)
    mean <- mean(QC$mu_per_read)
    Q30 <- (sum(QC$mu_per_read>=30)/length(QC$mu_per_read))*100

    scores_qc[i, "Q_Mean_R2"] <- mean
    scores_qc[i, "%_>=Q30_R2"] <- Q30

    i= i+1
  }
  #saveRDS(scores_qc, file= sprintf("%s/scores_QC.rds", patients_dir))
  trimeado <- ifelse(trimmed == TRUE, "trimmed", "")
  write.xlsx(scores_qc, file= paste(patients_dir, "/MetricasQC_", trimeado,".xlsx", sep =""), rowNames= TRUE)

  return(scores_qc)

}


patients_dir <- "/media/8tb02/Daniela/Leucemia"
bbmap = TRUE
mycofree = "YES"
path_bbmap = "/media/8tb02/Daniela/Leucemia/bbmap"
index_rsubread = "/media/8tb02/Daniela/Leucemia/Mycoplasma/indexRsubread/mycoplasma_index"
rate = "0.6"

library(dplyr)

# Contenido del archivo Rmd
rmd_content <- sprintf('
---
title: "Reporte de Calidad de Cohorte"
output:
  pdf_document: default
  html_document:
    df_print: paged
editor_options:
  markdown:
    wrap: 72
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(readr)
library(stringr)


#Métricas de Alineamiento ----------------------------------------
  #metricasSTAR <- calculateQCMetricsSTAR(patients_dir = "%s", trimmed = TRUE)
  metrics_long <- calculateQCMetricsSTAR(patients_dir = "%s", trimmed = TRUE)
  #ids_destacados <- c("21", "48", "15")



  metrics_long <- pivot_longer(
    metrics_long,
    cols = -ID,
    names_to = "Métrica",
    values_to = "Valor"
  )
  metrics_long <- mutate(
    metrics_long,
    ID = as.character(ID),
    etiqueta = ID,
    color_id = ID
  )

  # Colores únicos para cada ID destacado
  #colores_personalizados <- c(
  #  "21" = "blue",
  #  "48" = "#E41A1C",
  #  "15" = "#E41A1C",
  #  "Otros" = "grey70"
  #)
  #library(dplyr)

  metrics_long_filtrado <- group_by(metrics_long, Métrica)
  metrics_long_filtrado <- mutate(
    metrics_long_filtrado,
    Q1 = quantile(Valor, 0.25),
    Q3 = quantile(Valor, 0.75),
    IQR = Q3 - Q1,
    es_outlier = Valor < Q1 | Valor > Q3
  )



  alignment_plot <- ggplot(metrics_long, aes(x = "", y = Valor)) +
    geom_boxplot() +
    geom_violin(fill = "lightblue", width = 1, alpha = 0.4, color = NA) +
    #geom_jitter(aes(color = color_id), width = 0, alpha = 0.7, size = 2) +
    geom_jitter(color = "black", width = 0, alpha = 0.7, size = 1.5) +
    #geom_text(aes(label = etiqueta), hjust = -0.5, size = 3, na.rm = TRUE) +
    geom_text(
    data = subset(metrics_long_filtrado, es_outlier),
    aes(label = etiqueta),
    hjust = -0.5,
    size = 2,
    na.rm = TRUE
  ) +
    facet_wrap(~ Métrica, scales = "free_y", nrow = 2) +
    #scale_color_manual(values = colores_personalizados) +
    theme_minimal() +
    theme(
      #strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Métricas Alineamiento", y = "Valor", x = NULL)

  print(alignment_plot)

  # Metricas de secuenciación -------------------------------------



  library(qckitfastq)
  metricsSeq <- calculateQCMetricsSeq(patients_dir = "%s", trimmed = FALSE)
  colnames(metricsSeq)[which(colnames(metricsSeq) == "Sample")] <- "ID"
  ReporteFinal_QC <- merge(metricsSeq, metricasSTAR, by ="ID")

  seq_plot <- plotFastQC_PBSQ(patients_dir = "%s", trimmed = FALSE)
  print(seq_plot)

# Mycofree ---------------------------------------------------


if("%s" != "-") {
     mapeoMycoplasma_Grupal <- mapeoMycoplasma_Grupal(dir_patients = "%s",
                                                   bbmap = TRUE,
                                                   path_bbmap = "%s",
                                                   index_rsubread = "%s",
                                                   rate = "%s")
    library(ggplot2)

    myco_plot <- ggplot(mapeoMycoplasma_Grupal, aes(x = ID, y = Mapeo)) +
      geom_col(fill = "steelblue") +
      geom_hline(yintercept = 0, color = "green", size = 0.8) +  # línea sobre eje 0
      scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # límite mínimo 0
      theme_minimal() +
      labs(title = "Mapeo a Mycoplasma por muestra", y = "Cantidad de lecturas mapeadas", x = "ID") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    print(myco_plot)
    ReporteFinal_QC <- merge(ReporteFinal_QC, mapeoMycoplasma_Grupal, by ="ID")

} else {
    message("No se realiza el control de contaminación por mycoplasma.")
    myco_plot <- ggplot()
}


titulo_seq = "Métricas de Calidad de Sequenciación"
texto_seq = "Cada curva representa la calidad promedio de cada muestra. Para obtener una buena calidad de secuenciación, todas las curvas deberían estar por encima de la línea horizontal verde, es decir que su calidad debe ser >Q28."
titulo_align = "Métricas de Calidad de Alineamiento"
texto_align = paste0("Cada gráfico representa la distribución de los valores de métricas relevantes de alineamiento.", "\n\n",
"1. Promedio de largo de lectura de entrada: Depende del tamaño de la fragmentación de la sequenciación.", "\n\n",
"2. Promedio de largo de lecturas mapeadas: Debe ser similar al gráfico 1.",  "\n\n",
"3. Número de lecturas: ",  "\n\n",
"4. Porcentaje de mapeo único: Porcentaje de la cantidad de lecturas que mapearon al genoma humano. Se espera que sean >80.",  "\n\n",
"5. Número de lecturas de mapeo único: cantidad de lecturas que mapearon al genoma humano."
)

titulo_myco = "Control de contaminación con Mycoplasma"
texto_myco = paste0("En los cultivos celulares, la contaminación por Mycoplasma es un problema frecuente y persistente, ya que estas bacterias son muy pequeñas, carecen de pared celular y no generan turbidez visible en el medio, lo que dificulta su detección. Esta ausencia de pared celular las hace resistentes a antibióticos que actúan sobre su síntesis, como la penicilina, y aunque existen antibióticos efectivos, su uso continuo no se recomienda por posible citotoxicidad. La presencia de micoplasma puede alterar el metabolismo celular, modificar la expresión génica y comprometer la validez de los datos experimentales.", "\n\n",
"Por eso, en nuestro control de calidad de muestras de RNAseq provenientes de cultivos, incluimos la detección de micoplasma para garantizar la fiabilidad de los resultados.", "\n\n",
"Ref: https://doi.org/10.1093/nar/gkv136 ")

```

```{r, echo=FALSE}
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0("# ", titulo_seq, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(texto_seq, "\n\n"))
cat("\n\n")
print(seq_plot)

cat("\n\n")
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0("# ", titulo_align, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(texto_align, "\n\n"))
cat("\n\n")
print(alignment_plot)

cat("\n\n")
cat("\n\n")
cat("\n\n")
knitr::asis_output(paste0("# ", titulo_myco, "\n\n"))
cat("\n\n")
knitr::asis_output(paste0(texto_myco, "\n\n"))
cat("\n\n")
print(myco_plot)

```
  ', patients_dir, patients_dir, patients_dir, patients_dir, mycofree, patients_dir, path_bbmap, index_rsubread, rate)


