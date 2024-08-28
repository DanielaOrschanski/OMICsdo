#' @title Run FastQC
#' @description Executes the FastQC for R1 and R2.
#' @param patent_dir Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder.
#' @import stringr
#' @import viridis
#' @import reshape2
#' @import readr
#' @export
runFastQC <- function(patient_dir) {

  file_list <- list.files(patient_dir)

  #Para que se pueda poner como entrada la carpeta de los trimmeados o la carpeta original:
  if  (startsWith(basename(patient_dir), "trimmed")) {
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "val_1.fq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("val_1.fq%s", gzip))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("val_2.fq%s", gzip))], sep="")
  } else {
    gzip <- ifelse(length(nchar(file_list[endsWith(file_list, "R1.fastq.gz")])) == 0, "", ".gz")
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R1.fastq%s", gzip))], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, sprintf("R2.fastq%s", gzip))], sep="")
  }

  if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
    stop("There are no fastq files in this directory")
  }

  #Evita repetir el analisis si ya fue hecho
  if (startsWith(basename(patient_dir), "trimmed")) {
    if (!(length(nchar(file_list[endsWith(file_list, "val_1_fastqc")])) == 0) & !(length(nchar(file_list[endsWith(file_list, "val_2_fastqc")])) == 0)) {
      message("The FastQC for this sample has already been done.")
      return(paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1_fastqc")], sep=""))
    }
  } else {
    if ((length(nchar(file_list[endsWith(file_list, "R1_fastqc")])) != 0) & (length(nchar(file_list[endsWith(file_list, "R2_fastqc")])) != 0)) {
      message("The FastQC for this sample has already been done.")
      return(paste0(patient_dir, "/", file_list[endsWith(file_list, "R1_fastqc")], sep=""))
    }
  }

  #Ejecuta el FastQC R1
  system2(FastQC, fileR1)
  file_list <- list.files(patient_dir)

  #Extrae el zip R1
  file_fastqc_zip <- paste0(patient_dir, "/", file_list[endsWith(file_list, "1_fastqc.zip")], sep="")
  unzip(file_fastqc_zip, exdir = patient_dir)

  #Ejecuta el FastQC R2
  system2(FastQC, fileR2)
  file_list <- list.files(patient_dir)

  #Extrae el zip R2
  file_fastqc_zip <- paste0(patient_dir, "/", file_list[endsWith(file_list, "2_fastqc.zip")], sep="")
  unzip(file_fastqc_zip, exdir = patient_dir)

  message("FastQC's analysis has finished!")

}

#' @title plotFastQC_PBSQ
#' @description generates the principal plot ("Per Base Sequence Quality") which compares the quality for all the samples.
#' @param patients_dir Path of the directory that contains one folder with fastq files (R1 and R2) of each patient.
#' @param trimmed it is set to TRUE when the fastQC files that you want to plot came from a trimmed file.
#' @param R will indicate if the curve of the plot will be constructed by the mean of R1 (R= "R1"), R2 (R= "R2") or by both (R= "R1R2").
#' @export
#' @import readr
#' @import stringr
#' @import ggplot2
#' @import reshape2
#' @import png
plotFastQC_PBSQ <- function(patients_dir, trimmed = FALSE, R= "R1R2") {

  #Separación de los modulos del fastqc para generar graficos
  dir_list <- list.dirs(path = patients_dir, full.names = FALSE, recursive = FALSE)
  cant_patients <- length(dir_list)

  list_modulo_PBSQ <- list() #Va a almacenar los dataframes para el grafico Per Base Sequence Quality

  #Recorro cada paciente para guardar las tablas para el grafico y sacar el promedio de las medias.
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
  plot2 <- plot +
    geom_line(data = df.long, aes(x = Base, y = value, color = variable, group = variable), size = 1) +
    labs(title = "Per Base Sequence Quality - 1 curva por paciente", x = "Base", y = "Quality Value")

  indicador <- ifelse(trimmed == TRUE, "_trimmeado", "_SINtrimmeado")
  png(filename = sprintf("%s/plotFastQC_PBSQ%s.png", patients_dir, indicador), width = 800, height = 600)  # Ajusta el tamaño según tus necesidades
  print(plot2)
  dev.off()
  return(sprintf("The plot from FastQC was saved in %s/plotFastQC_PerBaseSequenceQuality.png", patients_dir))

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

