#' @title Obtain statistics from fusions.
#' @description Generates 2 dataframes: 1. Rbind from every fusion report of the samples and 2. Stats from all samples indicating group.
#' @param patients_dir path to the folder that contains one folder for each sample.
#' @param Metadata dataframe that contains information about all the samples. One sample per row.
#' @param group must be a name of one of the columns of the metadata
#' @return Todos_FusionReport and Stats_Fusions: 2 dataframes.
#' @export
#' @import openxlsx
#' @import readxl
#' @examples fusionStats(patients_dir, Metadata, group = "group")

fusionStats <- function(patients_dir, Metadata = NA, group = NA, cohorte = "") {

  ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  length(ids)
  #ids <- ids[-1]

  #MetadataSRA <- read.table("/media/4tb1/Daniela/Environ/Fusiones/SraRunTable.txt", header = TRUE, sep = ",")
  #colnames(MetadataSRA)[which(colnames(MetadataSRA) == "metastasis")] <- "MTT"


  Todos_FusionReport <- data.frame()
  Stats_Fusions <- data.frame()
  k=1
  #id <- ids[1]
  for (id in ids) {
    i <- basename(id)
    print(i)
    fusions_file <- sprintf("%s/trimmed/%s_FusionReport.xlsx", id, i)
    FusionReport <- read_excel(fusions_file)

    if(!is.na(Metadata)) {
      Grupo <- as.character(Metadata[which(Metadata$ID == i), group])
      FusionReport$MTT <- as.character(Metadata[which(Metadata$ID == i), "MTT"])
      FusionReport$Grupo <- Grupo

    }

    FusionReport$ID <-i
    FusionReport_ID <- cbind(FusionReport$ID, FusionReport[,1:(ncol(FusionReport)-1)])
    colnames(FusionReport)[1] <- "ID"
    FusionReport$Cohorte <- cohorte

    if(nrow(Todos_FusionReport) == 0) {
      Todos_FusionReport <- FusionReport_ID
    } else {
      Todos_FusionReport <- rbind(Todos_FusionReport, FusionReport_ID)
    }


    #Estadisticas de los reportes -----------------------------------

    Stats_Fusions[k, "ID"] <- i
    Stats_Fusions[k, "Cantidad_Fusiones"] <- nrow(FusionReport)

    if(!(nrow(FusionReport) == 0)) {
      contador_high <- 0
      contador_medium <- 0
      contador_low <- 0
      for (j in 1:nrow(FusionReport)) {
        if(FusionReport$confidence[j] == "high"){
          contador_high <- contador_high + 1
        } else if (FusionReport$confidence[j] == "medium"){
          contador_medium <- contador_medium + 1
        } else {
          contador_low <- contador_low + 1
        }
      }

      Stats_Fusions[k, "Fusiones_conf_H"] <- contador_high
      Stats_Fusions[k, "Fusiones_conf_M"] <- contador_medium
      Stats_Fusions[k, "Fusiones_conf_L"] <- contador_low
    } else {
      print(sprintf("Se encontraron 0 fusiones en %s", i))
    }


    if(!is.na(Metadata)) {
      met <- as.character(Metadata[which(Metadata$ID == i), "MTT"])
      Grupo <- as.character(Metadata[which(Metadata$ID == i), group])
      Stats_Fusions[k, "Grupo"] <- Grupo
      Stats_Fusions[k, "MTT"] <- met
    }
    Stats_Fusions$Cohorte <- cohorte

    k = k+1
  }


  Stats_Fusions[is.na(Stats_Fusions)] <- 0
  openxlsx::write.xlsx(as.data.frame(Todos_FusionReport), file = sprintf("%s/Todos-FusionReports_%s.xlsx", patients_dir, cohorte))
  write.xlsx(Stats_Fusions, file = sprintf("%s/StatsFusions_%s.xlsx", patients_dir, cohorte))

  TFB <- Todos_FusionReport[which(Todos_FusionReport$confidence == "high"), c(1,2,3,6,7, 16, 38, 37)]

  # Unir gene1 y gene2 en una sola columna y mantener información de la muestra
  TFB_long <- TFB %>%
    pivot_longer(cols = c(gene1, gene2), names_to = "Gene_Type", values_to = "Gene")

  str(TFB_long)
  colnames(TFB_long)[1] <- "ID"
  TFB_long$ID <- factor(TFB_long$ID)

  # Contabilizar métricas por gen
  if(!is.na(group)){
    categorias_grupo <- unique(Stats_Fusions$Grupo)

    gen_counts <- TFB_long %>%
      group_by(Gene) %>%
      summarise(
        Total_Apariciones = n(),  # Cantidad total de veces que aparece el gen
        Muestras_Distintas = n_distinct(ID),  # Muestras únicas donde aparece
        MET_Pos = sum(MTT == "MET+"),  # Veces que aparece en muestras MET+,
        MET_Neg = sum(MTT == "MET-"),
        LumHER2 = sum(Grupo == "Luminal HER2"),
        LumB = sum(Grupo == "Luminal B"),
        TripleNeg = sum(Grupo == "Triple Negativo"),
        HER2 = sum(Grupo == "HER2"),
        LumA = sum(Grupo == "Luminal A")
        #Apariciones_TumorTissue = sum(Grupo == "breast tumor")
      ) %>%
      arrange(desc(Total_Apariciones))  # Ordenar por frecuencia

  }

  gen_counts <- TFB_long %>%
    group_by(Gene) %>%
    summarise(
      Total_Apariciones = n(),  # Cantidad total de veces que aparece el gen
      Muestras_Distintas = n_distinct(ID),  # Muestras únicas donde aparece
      Apariciones_MET_Pos = sum(MTT == "MET+"),  # Veces que aparece en muestras MET+,
      Apariciones_MET_Neg = sum(MTT == "MET-")
      #,  # Veces que aparece en muestras MET+
      #Apariciones_NormalTissue = sum(Grupo == "normal breast tissue"),
      #Apariciones_TumorTissue = sum(Grupo == "breast tumor")
    ) %>%
    arrange(desc(Total_Apariciones))  # Ordenar por frecuencia

  write.xlsx(gen_counts, file = sprintf("%s/GeneFusions_%s.xlsx", patients_dir, cohorte))

  #Generar boxplot:
  boxplots_TFB_MTT(stats = Stats_Fusions, group = group, cohorte = cohorte)

  #Generar analisis sobrevida:
  analisis_sobrevida(stats = Stats_Fusions, metadata = Metadata)


  return(list(Todos_FusionReport, Stats_Fusions, gen_counts))

}


boxplots_TFB_MTT <- function(stats, group, cohorte) {

  max <- max(stats$Fusiones_conf_H)
  paso <- round(max/10)

  box_TFB_MTT <- ggplot(stats, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
    labs(title = sprintf("TFB por MTT - %s", cohorte),
         x = "MTT",
         y = "Cantidad de Fusiones") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, max +paso, by = paso)) +
    stat_compare_means(method = "wilcox.test", label = "p.format")  # Agrega p-valor con Wilcoxon

  print(box_TFB_MTT)

  if(!is.na(group)) { #Boxplot por subtipo:
    box_TFB_MTT_Grupo <- ggplot(stats, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
      labs(title = sprintf("TFB por MTT y Grupo - %s", cohorte),
           x = "MTT",
           y = "Cantidad de Fusiones") +
      theme_minimal() +
      #scale_y_continuous(limits = c(0, NA)) +
      scale_y_continuous(limits = c(0, NA), breaks = seq(0, max+paso, by = paso)) +
      stat_compare_means(method = "wilcox.test", label = "p.format") + # Agrega p-valor con Wilcoxon
      facet_grid( ~ Grupo)  # Facet por Cohorte y Subtipo

    print(box_TFB_MTT_Grupo)
  }
}

#Necesito que metadata tenga info de tiempo y MTT:

analisis_sobrevida <- function(stats, metadata) {

  thr <- median(stats$Fusiones_conf_H)
  #thr <- mean(stats$Fusiones_conf_H)

  stats$group <- ifelse(stats$Fusiones_conf_H > thr, "High", "Low")
  stats$group <- factor(stats$group, levels = c("Low","High"))

  #incorporar informacion de sobrevida o sobrevida libre de progresion:
  stats_con_tiempo <- merge(stats, metadata[,c("ID", "tiempo")], by = "ID", all = FALSE)
  stats_con_tiempo$tiempo <- as.numeric(stats_con_tiempo$tiempo)
  if(any(is.na(stats_con_tiempo$tiempo))) {
    stats_con_tiempo <- stats_con_tiempo[!is.na(stats_con_tiempo$tiempo),]
  }
  stats_con_tiempo$evento <- ifelse(stats_con_tiempo$MTT == "MET-", 0, 1)
  library("survival")
  library("survminer")

  surv_object <- Surv(time = stats_con_tiempo$tiempo, event = stats_con_tiempo$evento)
  fit1 <- survfit(surv_object ~ group, data = stats_con_tiempo)

  curva_sobrevida <- ggsurvplot(fit1, data = stats_con_tiempo, size = 1,  # change line size
                                linetype = "strata", # change line type by groups
                                palette = c("red","blue"), # custom color palette
                                conf.int = TRUE, # Add confidence interval
                                pval = TRUE, # Add p-value
                                risk.table = TRUE, # add table
                                title = sprintf("Curva Sobrevida Libre de Progresión - %s - thr = %s", cohorte, thr)
  )

  print(curva_sobrevida)
  table_matrix <- table(stats_con_tiempo$group, stats_con_tiempo$MTT)
  print(table_matrix)
}
