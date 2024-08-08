#' @title Calculate QC metrics - STAR
#' @description Executes tkjvnfkd.
#' @param patients_dir Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder.
#' @param trimmed kjvdfk
#' @return metricasSTAR df with metrics from alignment
#' @export
calculateQCMetricsSTAR <- function(patients_dir, trimmed = FALSE) {

  dir_list <- list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  dir_list <- dir_list[-1]
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

  write.xlsx(metricasST, file = sprintf("%s/MetricasQCSTAR.xlsx", patients_dir))
  return(metricasSTAR)

}
