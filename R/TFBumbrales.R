#Probar dif umbrales de TFB en environ y dar % de MTT para cada uno
# varuar el umbral desde 0 a maximo de a 1
# indep de subtipo y para cada subtipo

Environ_Paula_metadata_global_2025 <- read_excel("/media/16TBDisk/Daniela/Environ/Environ_Paula_metadata_global_2025.xlsx")

stats_CDG <- as.data.frame(read_excel("/media/4tb1/Paula/13Br_Validacion/Muestras/StatsFusions.xlsx"))
stats_tcl <- as.data.frame(read_excel("/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/Muestras/StatsFusions.xlsx"))
library(openxlsx)
stats_mama61 <- as.data.frame(read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/StatsFusions.xlsx"))

#metadata_mama61 <- Environ_Paula_metadata_global_2025[which(Environ_Paula_metadata_global_2025$categoria_1 == "mama_61"),]
#colnames(metadata_mama61)
#metadata_mama61$ID[which(metadata_mama61$ID == "Br0310")] <- "Br0310-1"
#metadata_mama61$ID[which(metadata_mama61$ID == "Br0307")] <- "Br0307-2"
#out_mama61 <- fusionStats(patients_dir = "/media/16TBDisk/Daniela/Environ/TodoMama", Metadata = metadata_mama61, group = "subtipo")
#stats_mama61 <- out_mama61[[2]]

#TCGA ---------------------------------------------------------------
#library(qtl2)
BRCA_TCGA <- read.csv("/media/4tb1/Paula/BRCA_TCGA.csv")
colnames(BRCA_TCGA)[2] <- "ID"
colnames(BRCA_TCGA)[3] <- "Fusiones_conf_H"
colnames(BRCA_TCGA)[which(colnames(BRCA_TCGA) == "ajcc_pathologic_m")] <- "MTT"
BRCA_TCGA <- BRCA_TCGA[which(BRCA_TCGA$MTT %in% c("M0", "M1")),]
table(BRCA_TCGA$MTT)
BRCA_TCGA$MTT <- ifelse(BRCA_TCGA$MTT == "M0", "MET-", "MET+")
stats_TCGA <- BRCA_TCGA[, c("ID", "Fusiones_conf_H", "MTT")]
str(stats_TCGA)

metadata_tcga <- readRDS("/media/4tb2/Daniela/Environ/FusionesActualizado/BRCA_rna.rds")
metadata_tcga <- metadata_tcga$targets
colnames(stats_TCGA)
library(dplyr)
stats_TCGA <- stats_TCGA %>%
  filter(MTT == "MET+" | (MTT == "MET-" & Grupo %in% c("Her2", "LumB")))

ggplot(stats_TCGA, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "Distribución de la Cantidad de Fusiones HIGH por MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))

# -------------------------------------------------
stats_mama61$Cohorte <- "Mama61"
stats_CDG$Cohorte <- "CDGenomics"
stats_tcl$Cohorte <- "TCL"
stats_TCGA$Cohorte <- "TCGA"
#stats_52SRA$Cohorte <- "SRA"
colnames(stats_mama61)
stats_conjunto <- rbind(stats_mama61[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_CDG[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_tcl[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        #stats_52SRA[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_TCGA)


TFB_umbrales_mama61 <- generar_tabla_TFB_umbrales(DB = stats_mama61)
plot_mama61 <- generar_TFB_umbrales_porSubtipo(DB = stats_mama61)
tabla_mama61 <- TFB_umbrales_mama61[[1]]
colnames(tabla_mama61)

TFB_umbrales_CDG <- generar_tabla_TFB_umbrales(DB = stats_CDG)

TFB_umbrales_tcl <- generar_tabla_TFB_umbrales(DB = stats_tcl)
plot_tcl <- generar_TFB_umbrales_porSubtipo(DB = stats_tcl)

TFB_umbrales_TCGA <- generar_tabla_TFB_umbrales(DB = stats_TCGA)
plot_tcga <- generar_TFB_umbrales_porSubtipo(DB = stats_TCGA)

wb <- createWorkbook()
addWorksheet(wb, "Mama61")
writeData(wb, "Mama61", TFB_umbrales_mama61)
addWorksheet(wb, "CDGenomics")
writeData(wb, "CDGenomics", TFB_umbrales_CDG)
addWorksheet(wb, "TCL")
writeData(wb, "TCL", TFB_umbrales_tcl)
addWorksheet(wb, "TCGA")
writeData(wb, "TCGA", TFB_umbrales_TCGA)

saveWorkbook(wb, "/media/4tb2/Daniela/Environ/FusionesActualizado/Fusiones-paraEnviron-2025/TFB_umbrales_indep.xlsx", overwrite = TRUE)



#-----------------------------------------------



generar_TFB_umbrales_porSubtipo <- function(DB) {

  subtipos <- unique(DB$Grupo)
  #subtipo = subtipos[2]
  list_tablas <- list()
  list_curvas <- list()
  list_rocs <- list()

  for (subtipo in subtipos) {
    print(subtipo)
    DB_subset <- DB[DB$Grupo == subtipo,]
    out <- generar_tabla_TFB_umbrales(DB = DB_subset)
    TFB_umb <- out[[1]]
    curva <- out[[2]]
    roc <- out[[3]]
    list_tablas[[subtipo]] <- TFB_umb
    list_curvas[[subtipo]] <- curva
    list_rocs[[subtipo]] <- roc
  }

  wb <- createWorkbook()

  # Iterar sobre la lista y agregar cada dataframe en una hoja con su nombre
  for (nombre in names(list_tablas)) {
    addWorksheet(wb, nombre)  # Agregar una hoja con el nombre del subtipo
    writeData(wb, nombre, list_tablas[[nombre]])  # Escribir el dataframe
  }

  cohorte <- unique(DB$Cohorte)
  saveWorkbook(wb, sprintf("/media/4tb2/Daniela/Environ/FusionesActualizado/Fusiones-paraEnviron-2025/TFB_umbrales_porSubtipo_%s.xlsx", cohorte), overwrite = TRUE)

  # Reducir tamaño del título en cada gráfico
  list_curvas <- lapply(list_curvas, function(p) {
    p + theme(plot.title = element_text(size = 10))  # Ajusta el tamaño del título
  })
  list_rocs <- lapply(list_rocs, function(p) {
    p + theme(plot.title = element_text(size = 10))  # Ajusta el tamaño del título
  })

  library(gridExtra)
  grid <- grid.arrange(grobs = list_curvas, ncol = 2)
  print(grid)
  grid_roc <- grid.arrange(grobs = list_rocs, ncol = 2)
  print(grid_roc)

  return(list(list_tablas, grid, grid_roc))

}

DB = stats_TCGA

generar_tabla_TFB_umbrales <- function(DB) {

  colnames(DB)[colnames(DB) == "Fusiones_conf_H"] <- "TFB"

  res <- data.frame()
  for (i in 0:max(DB$TFB)) {
    thr <- i  # Set threshold
    DB$group <- ifelse(DB$TFB > thr, "TFB_H", "TFB_L")  # Classify based on threshold

    #Para hacer lo del score:
    #DB$Score <- ifelse(DB$TFB <= thr, "TFB_H", "TFB_L")

    # Count occurrences for each category
    count_TFB_H_Pos <- sum(DB$group == "TFB_H" & DB$MTT == "MET+")
    count_TFB_H_Neg <- sum(DB$group == "TFB_H" & DB$MTT == "MET-")
    count_TFB_L_Pos <- sum(DB$group == "TFB_L" & DB$MTT == "MET+")
    count_TFB_L_Neg <- sum(DB$group == "TFB_L" & DB$MTT == "MET-")

    # Calculate total count
    total_count <- count_TFB_H_Pos + count_TFB_H_Neg + count_TFB_L_Pos + count_TFB_L_Neg

    # Calculate percentages
    percent_TFB_H_Pos <- ifelse(total_count > 0, (count_TFB_H_Pos / total_count) * 100, 0)
    percent_TFB_H_Neg <- ifelse(total_count > 0, (count_TFB_H_Neg / total_count) * 100, 0)
    percent_TFB_L_Pos <- ifelse(total_count > 0, (count_TFB_L_Pos / total_count) * 100, 0)
    percent_TFB_L_Neg <- ifelse(total_count > 0, (count_TFB_L_Neg / total_count) * 100, 0)


    #Recalcular los % pero / cant_MTT
    count_pos <-  count_TFB_H_Pos + count_TFB_L_Pos
    percent_TFB_L_Pos_SCORE <- ifelse(count_pos > 0, (count_TFB_L_Pos / count_pos) * 100, 0)


    # Append to results
    res <- rbind(res, data.frame(
      thr,
      count_TFB_H_Pos, count_TFB_H_Neg, count_TFB_L_Pos, count_TFB_L_Neg,
      percent_TFB_H_Pos, percent_TFB_H_Neg, percent_TFB_L_Pos, percent_TFB_L_Neg,
      percent_TFB_L_Pos_SCORE
    ))
  }

  #Graficas cada TFB_umbral:
  TFB_umbrales <- res

  cohorte <- unique(DB$Cohorte)

  #Le pongo el nombre de subtipo si es que el analisis se esta haciendo por subtipo, sino pongo indep
  subtipo <- unique(DB$Grupo)
  if(length(subtipo) >1 | length(subtipo) == 0) {
    subtipo = "indep"
  }

  #Armar una tabla tal que para cada umbral haya un par x(FPR), y(TPR)
  # Es como si pensaramos una matriz de confusion donde:
  #       TSB_H     TSB_L
  # MTT +
  # MTT -

  roc_data <- TFB_umbrales %>%
    mutate(
      TPR = count_TFB_H_Pos / (count_TFB_H_Pos + count_TFB_L_Pos),
      FPR = count_TFB_H_Neg / (count_TFB_H_Neg + count_TFB_L_Neg)
    )

  # Graficar la curva ROC
  library(ggrepel)
  roc <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "red") +
    geom_point() +
    geom_text_repel(aes(label = thr), size = 3, box.padding = 0.3, max.overlaps = 20) +  # Usar geom_text_repel

    #geom_text(aes(label = thr), vjust = 1.5, hjust = 1, size = 3) +  # Agregar el valor del umbral
    labs(title = sprintf("Curva ROC - %s - %s", cohorte, subtipo),
         x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)") +
    theme_minimal()+
    theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))  # Aumentar margen superior si es necesario

  print(roc)

  library(ggplot2)
  library(tidyr)
  TFB_umbrales_long <- pivot_longer(TFB_umbrales, cols = starts_with("percent"),
                                    names_to = "Categoria", values_to = "Valor")

  colores_personalizados <- c(
    "percent_TFB_H_Pos" = "#E41A1C",  # Rojo fuerte
    "percent_TFB_H_Neg" = "#FF5733",  # Rojo anaranjado
    "percent_TFB_L_Pos" = "#377EB8",  # Azul fuerte
    "percent_TFB_L_Neg" = "#6A9CE6"   # Azul más claro
  )
  curva <- ggplot(TFB_umbrales_long, aes(x = thr, y = Valor, color = Categoria)) +
    geom_point() +  # Agrega puntos
    geom_line() +   # Agrega líneas entre los puntos
    scale_color_manual(values = colores_personalizados) +
    labs(title = sprintf("Porcentaje de TFB según umbral - %s - %s", cohorte, subtipo),
         x = "Threshold (thr)",
         y = "Porcentaje",
         color = "Categoría") +
    theme_minimal()
  print(curva)

  return(list(res, curva, roc))
}
