#' @title Obtain statistics from fusions.
#' @description Generates 2 dataframes: 1. Rbind from every fusion report of the samples and 2. Stats from all samples indicating group.
#' @param patients_dir path to the folder that contains one folder for each sample.
#' @param Metadata dataframe that contains information about all the samples. One sample per row.
#' @param group must be a name of one of the columns of the metadata
#' @return Todos_FusionReport and Stats_Fusions: 2 dataframes.
#' @export
#' @import openxlsx
#' @import readxl
#' @import ggpubr
#' @examples fusionStats(patients_dir, Metadata, group = "group")

fusionStats <- function(patients_dir, Metadata = NA, group = NA, cohorte = "", sobrevida = TRUE, Apareados = FALSE, grupo_ids_apareados = NA) {

  ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
  length(ids)
  #ids <- ids[-1]
  if(any(grepl(".git", ids))) {
    ids <- ids[-which(grepl(".git", ids))]
  }

  #MetadataSRA <- read.table("/media/4tb1/Daniela/Environ/Fusiones/SraRunTable.txt", header = TRUE, sep = ",")
  #colnames(MetadataSRA)[which(colnames(MetadataSRA) == "metastasis")] <- "MTT"

  if(!all(is.na(Metadata))) {
    ids_analizar <- Metadata$ID
    ids <- ids[which(basename(ids) %in% ids_analizar)]
  }


  Todos_FusionReport <- data.frame()
  Stats_Fusions <- data.frame()
  k=1
  #id <- ids[1]
  for (id in ids) {
    i <- basename(id)
    print(i)
    fusions_file <- sprintf("%s/trimmed/%s_FusionReport.xlsx", id, i)
    FusionReport <- read_excel(fusions_file)

    #if(!missing(Metadata)) {
    #  Grupo <- as.character(Metadata[which(Metadata$ID == i), group])
    #  FusionReport$MTT <- as.character(Metadata[which(Metadata$ID == i), "MTT"])
    #  FusionReport$Grupo <- Grupo
    #}

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


    #if(!missing(Metadata)) {
    #  met <- as.character(Metadata[which(Metadata$ID == i), "MTT"])
    #  Grupo <- as.character(Metadata[which(Metadata$ID == i), group])
    #  Stats_Fusions[k, "Grupo"] <- Grupo
    #  Stats_Fusions[k, "MTT"] <- met
    #}
    Stats_Fusions$Cohorte <- cohorte

    k = k+1
  }


  if(any(is.na(Stats_Fusions))) {
    Stats_Fusions[is.na(Stats_Fusions)] <- 0
  }

  if(!all(is.na(Metadata))) {
    Stats_Fusions <- merge(Stats_Fusions, Metadata, by = "ID")
  }

  openxlsx::write.xlsx(as.data.frame(Todos_FusionReport), file = sprintf("%s/Todos-FusionReports_%s.xlsx", patients_dir, cohorte))

  # Unir gene1 y gene2 en una sola columna y mantener información de la muestra
  TFB_long <- Todos_FusionReport %>%
    pivot_longer(cols = c(gene1, gene2), names_to = "Gene_Type", values_to = "Gene")

  str(TFB_long)
  TFB_long <- as.data.frame(TFB_long)
  colnames(TFB_long)[1] <- "ID"
  TFB_long$ID <- factor(TFB_long$ID)

  # Contabilizar métricas por gen
  library(dplyr)
  gen_counts <- TFB_long %>%
    group_by(Gene) %>%
    summarise(
      Total_Apariciones = n(),  # Cantidad total de veces que aparece el gen
      Muestras_Distintas = n_distinct(ID)  # Muestras únicas donde aparece
    ) %>%
    arrange(desc(Total_Apariciones))

  if(!is.na(group) & !is.null(Stats_Fusions$MTT)) {
    categorias_grupo <- unique(Stats_Fusions$Grupo)
    if(group == "tissue_type") {
      gen_counts <- TFB_long %>%
        group_by(Gene) %>%
        summarise(
          Total_Apariciones = n(),  # Cantidad total de veces que aparece el gen
          Muestras_Distintas = n_distinct(ID),  # Muestras únicas donde aparece
          MET_Pos = n_distinct(ID[MTT == "MET+"]),
          MET_Neg = n_distinct(ID[MTT == "MET-"]),
          Adyacente = n_distinct(ID[Grupo == "Normal appearing thyroid tissue adjacent to PTC tumor"]),
          Normal = n_distinct(ID[Grupo == "Normal"]),
          Tumoral = n_distinct(ID[Grupo == "PTC tumor tissue"])
        ) %>%
        arrange(desc(Total_Apariciones))

    } else if(group == "Subtipo" | group == "subtipo") {
      gen_counts <- TFB_long %>%
        group_by(Gene) %>%
        summarise(
          Total_Apariciones = n(),  # Cantidad total de veces que aparece el gen
          Muestras_Distintas = n_distinct(ID),  # Muestras únicas donde aparece
          MET_Pos = n_distinct(ID[MTT == "MET+"]),  # Veces que aparece en muestras MET+,
          MET_Neg = n_distinct(ID[MTT == "MET-"]),
          LumHER2 = n_distinct(ID[Grupo == "Luminal HER2"]),
          LumB = n_distinct(ID[Grupo == "Luminal B"]),
          TripleNeg = n_distinct(ID[Grupo == "Triple Negativo"]),
          HER2 = n_distinct(ID[Grupo == "HER2"]),
          LumA = n_distinct(ID[Grupo == "Luminal A"])
        ) %>%
        arrange(desc(Total_Apariciones))  # Ordenar por frecuencia

    } else if(group == "Grupo" | group == "grupo") {
      gen_counts <- TFB_long %>%
        group_by(Gene) %>%
        summarise(
          Total_Apariciones = n(),  # Cantidad total de veces que aparece el gen
          Muestras_Distintas = n_distinct(ID),  # Muestras únicas donde aparece
          MET_Pos = n_distinct(ID[MTT == "MET+"]),  # Veces que aparece en muestras MET+,
          MET_Neg = n_distinct(ID[MTT == "MET-"]),
          GI = n_distinct(ID[Grupo == "Grupo I"]),
          GII = n_distinct(ID[Grupo == "Grupo II"]),
          GIII = n_distinct(ID[Grupo == "Grupo III"]),
          GIV = n_distinct(ID[Grupo == "Grupo IV"])
        ) %>%
        arrange(desc(Total_Apariciones))
    }
  }
  library(openxlsx)
  write.xlsx(gen_counts, file = sprintf("%s/GeneFusions_%s.xlsx", patients_dir, cohorte))

  #Generar boxplot:
  if(!all(is.na(Metadata))) {
    boxplots_TFB_MTT(stats = Stats_Fusions, group = group, cohorte = cohorte, Apareados = Apareados, grupo_ids_apareados = grupo_ids_apareados)
  }

  #Generar analisis sobrevida:
  if(sobrevida == TRUE) {
    analisis_sobrevida(stats = Stats_Fusions, metadata = Metadata, thre = thre)
  }

  #Calcular % de fusiones que son kinasas - Por paciente:
  library(dplyr)
  colnames(Todos_FusionReport)[1] <- "ID"

  total_kinase_fus <- Todos_FusionReport %>%
    group_by(ID) %>%
    summarise(total_kinase_fus = sum(grepl("Protein_kinase_domain", retained_protein_domains, ignore.case = TRUE))) %>%
    ungroup()

  high_kinase_fus <- Todos_FusionReport %>%
    filter(confidence == "high") %>%
    group_by(ID) %>%
    summarise(high_kinase_fus = sum(grepl("Protein_kinase_domain", retained_protein_domains, ignore.case = TRUE))) %>%
    ungroup()

  Stats_Fusions <- merge(Stats_Fusions, total_kinase_fus, by = "ID", all.x = TRUE)
  Stats_Fusions <- merge(Stats_Fusions, high_kinase_fus, by = "ID", all.x = TRUE)

  colnames(Stats_Fusions)
  Stats_Fusions$PercentageKinaseHigh <- (Stats_Fusions$high_kinase_fus / Stats_Fusions$Fusiones_conf_H )*100
  Stats_Fusions$PercentageKinase <- (Stats_Fusions$total_kinase_fus / Stats_Fusions$Cantidad_Fusiones )*100

  write.xlsx(Stats_Fusions, file = sprintf("%s/StatsFusions_%s.xlsx", patients_dir, cohorte))

  return(list(Todos_FusionReport, Stats_Fusions, gen_counts))

}


boxplots_TFB_MTT <- function(stats, group, cohorte, Apareados = FALSE, grupo_ids_apareados = NA) {

  library(ggpubr)
  max <- max(stats$Fusiones_conf_H)
  paso <- round(max/10)
  #stats$MTT <- as.factor(stats$MTT)

  while (!is.null(dev.list())) dev.off()

  if(!is.null(stats$MTT)) {

    # GRAFICO SEGUN MTT: -----------------------------------------------
    if(Apareados == FALSE) {
      test_result <- wilcox.test(Fusiones_conf_H ~ MTT, data = stats)
    } else if (Apareados == TRUE) {
      test_result <- wilcox.test(Fusiones_conf_H ~ MTT, data = stats, paired = TRUE)
    }

    p_value <- test_result$p.value
    p_label <- sprintf("p = %.3g", p_value)

    box_TFB_MTT <- ggplot(stats, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
      labs(title = sprintf("FPR por MTT - %s", cohorte),
           x = "MTT",
           y = "FPR") +
      theme_minimal() +
      scale_y_continuous(limits = c(0, NA), breaks = seq(0, max + paso, by = paso)) +
      annotate("text", x = 1,
               y = max(stats$Fusiones_conf_H, na.rm = TRUE) * 1.05,
               label = p_label, hjust = 0, size = 4)

    print(box_TFB_MTT)

    if(!is.na(group)) {
      # GRAFICO SEGUN MTT + GRUPO: -----------------------------------------------

      library(dplyr)

      stats$Grupo <- stats[[group]]
      #elimino los NA si hay:
      if(any(is.na(stats$Grupo))) {
        stats <- stats[-which(is.na(stats$Grupo)),]
      }
      stats$Grupo <- as.factor(stats$Grupo)

      wilcox_results <- stats %>%
        group_by(Grupo) %>%
        summarise(
          p_value = if (n_distinct(MTT) == 2) {
            wilcox.test(Fusiones_conf_H ~ MTT)$p.value
          } else {
            NA
          } )

      print(wilcox_results)

      box_TFB_MTT_Grupo <- ggplot(stats, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
        geom_violin(alpha = 0.5) +
        geom_boxplot() +
        geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
        labs(title = sprintf("FPR por MTT y Grupo - %s", cohorte),
             x = "MTT",
             y = "FPR") +
        theme_minimal() +
        scale_y_continuous(limits = c(0, NA), breaks = seq(0, max + paso, by = paso)) +
        facet_grid(~ Grupo)

      print(box_TFB_MTT_Grupo)
    }
  }

  if(!is.na(group)) {

    stats$Grupo <- stats[[group]]
    #elimino los NA si hay:
    if(any(is.na(stats$Grupo))) {
      stats <- stats[-which(is.na(stats$Grupo)),]
    }
    stats$Grupo <- as.factor(stats$Grupo)
    table(stats$Grupo)

    if( length(unique(stats$Grupo))>2 ) { #Kruskall
      kruskal_result <- kruskal.test(Fusiones_conf_H ~ Grupo, data = stats)
      p_value_FPR <- kruskal_result$p.value

      kruskal_result <- kruskal.test(Cantidad_Fusiones ~ Grupo, data = stats)
      p_value_Fus <- kruskal_result$p.value

    } else if ( length(unique(stats$Grupo)) == 2 ) { #Wilcox ------------

      if(Apareados == FALSE) {
        wilcox_result_FPR <- wilcox.test(Fusiones_conf_H ~ Grupo, data = stats)
        p_value_FPR <- wilcox_result_FPR$p.value

        wilcox_result_Fus <- wilcox.test(Cantidad_Fusiones ~ Grupo, data = stats)
        p_value_Fus <- wilcox_result_Fus$p.value

      } else if (Apareados == TRUE) {

        #APAREADOS - FPR
        stats$isolate <- stats[[grupo_ids_apareados]]
        stats_wide <- pivot_wider(
          stats,
          id_cols = isolate, #nombre que define los ID apareados
          names_from = Grupo,
          values_from = Fusiones_conf_H
        )
        wilcox_result_FPR <- wilcox.test(stats_wide[[2]], stats_wide[[3]], paired = TRUE)
        p_value_FPR <- wilcox_result_FPR$p.value

        #APAREADOS - Fus
        stats_wide <- pivot_wider(
          stats,
          id_cols = isolate, #nombre que define los ID apareados
          names_from = Grupo,
          values_from = Cantidad_Fusiones
        )
        wilcox_result_Fus <- wilcox.test(stats_wide[[2]], stats_wide[[3]], paired = TRUE)
        p_value_Fus <- wilcox_result_Fus$p.value
      }
      #wilcox_result <- wilcox.test(Fusiones_conf_H ~ Grupo, data = stats)
      #p_value <- wilcox_result$p.value
      #-------------------------------------------
    }

    p_label_FPR <- sprintf("p = %.3g", p_value_FPR)
    p_label_Fus <- sprintf("p = %.3g", p_value_Fus)

    box_TFB_MTT <- ggplot(stats, aes(x = Grupo, y = Fusiones_conf_H, fill = Grupo)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
      labs(title = sprintf("FPR por %s - %s", group, cohorte),
           x = sprintf("%s", group),
           y = "FPR") +
      theme_minimal() +
      #scale_y_continuous(limits = c(0, NA), breaks = seq(0, max + paso, by = paso)) +
      annotate("text", x = 1, y = max(stats$Fusiones_conf_H, na.rm = TRUE) * 1.05,
               label = p_label_FPR, hjust = 0, size = 4)

    print(box_TFB_MTT)

    box_Fus_MTT <- ggplot(stats, aes(x = Grupo, y = Cantidad_Fusiones, fill = Grupo)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot() +
      geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
      labs(title = sprintf("Fus por %s - %s", group, cohorte),
           x = sprintf("%s", group),
           y = "Fus") +
      theme_minimal() +
      #scale_y_continuous(limits = c(0, NA), breaks = seq(0, max + paso, by = paso)) +
      annotate("text", x = 1, y = max(stats$Fusiones_conf_H, na.rm = TRUE) * 1.05,
               label = p_label_Fus, hjust = 0, size = 4)

    print(box_Fus_MTT)

    if(Apareados == TRUE) {

      box_TFB_MTT <- ggplot(stats, aes(x = Grupo, y = Fusiones_conf_H, fill = Grupo)) +
        geom_violin(alpha = 0.5) +
        geom_boxplot() +
        geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
        # Añadir líneas punteadas entre puntos con el mismo isolate
        geom_line(aes(group = isolate), color = "gray10", linetype = "dotted", alpha = 0.6) +
        labs(title = sprintf("FPR por %s - %s", group, cohorte),
             x = sprintf("%s", group),
             y = "FPR") +
        theme_minimal() +
        annotate("text", x = 1, y = max(stats$Fusiones_conf_H, na.rm = TRUE) * 1.05,
                 label = p_label_FPR, hjust = 0, size = 4)

      print(box_TFB_MTT)

      box_Fus_MTT <- ggplot(stats, aes(x = Grupo, y = Cantidad_Fusiones, fill = Grupo)) +
        geom_violin(alpha = 0.5) +
        geom_boxplot() +
        geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
        # Añadir líneas punteadas entre puntos con el mismo isolate
        geom_line(aes(group = isolate), color = "gray10", linetype = "dotted", alpha = 0.6) +
        labs(title = sprintf("Fus por %s - %s", group, cohorte),
             x = sprintf("%s", group),
             y = "Fus") +
        theme_minimal() +
        annotate("text", x = 1, y = max(stats$Fusiones_conf_H, na.rm = TRUE) * 1.05,
                 label = p_label_Fus, hjust = 0, size = 4)

      print(box_Fus_MTT)

    }
  }

}

#Necesito que metadata tenga info de tiempo y MTT:
#stats <- Stats_Fusions
#metadata <- Metadata
analisis_sobrevida <- function(stats, metadata, thre = "median") {

  stats$Fusiones_conf_H <- as.numeric(stats$Fusiones_conf_H)

  thr <- ifelse(thre == "mean",
                mean(stats$Fusiones_conf_H),
                median(stats$Fusiones_conf_H))
  #thr <- mean(stats$Fusiones_conf_H)

  stats$group <- ifelse(stats$Fusiones_conf_H > thr, "High", "Low")
  stats$group <- factor(stats$group, levels = c("Low","High"))

  cohorte <- unique(stats$Cohorte)

  #incorporar informacion de sobrevida o sobrevida libre de progresion:
  #stats_con_tiempo <- merge(stats, metadata[,c("ID", "tiempo")], by = "ID", all = FALSE)
  stats_con_tiempo <- stats[, c("ID", "tiempo", "evento", "Fusiones_conf_H", "Cohorte", "group")]

  #Chequeo TIEMPO:
  stats_con_tiempo$tiempo <- as.numeric(stats_con_tiempo$tiempo)
  if(any(is.na(stats_con_tiempo$tiempo))) {
    stats_con_tiempo <- stats_con_tiempo[!is.na(stats_con_tiempo$tiempo),]
  }

  #Chequeo evento si es MTT o ya viene con evento:
  if(!is.null(stats_con_tiempo$MTT)) {
    stats_con_tiempo$evento <- ifelse(stats_con_tiempo$MTT == "MET-", 0, 1)
  }
  if(any(is.na(stats_con_tiempo$evento))) {
    stats_con_tiempo <- stats_con_tiempo[!is.na(stats_con_tiempo$evento),]
  }

  stats_con_tiempo$evento <- as.numeric(stats_con_tiempo$evento)

  library("survival")
  library("survminer")

  stats_con_tiempo$group <- factor(stats_con_tiempo$group, levels = c("High","Low"))

  surv_object <- Surv(time = stats_con_tiempo$tiempo, event = stats_con_tiempo$evento)
  fit1 <- survfit(surv_object ~ group, data = stats_con_tiempo)


  if(!is.null(stats$MTT)) {
    tipo_sobrevida <- "Libre de Metástasis"
  } else {
    tipo_sobrevida <- "OS"
  }
  curva_sobrevida <- ggsurvplot(fit1, data = stats_con_tiempo, size = 1,  # change line size
                                linetype = "strata", # change line type by groups
                                palette = c("red","blue"), # custom color palette
                                conf.int = TRUE, # Add confidence interval
                                pval = TRUE, # Add p-value
                                risk.table = TRUE, # add table
                                title = sprintf("Curva Sobrevida %s - %s - thr = %s", tipo_sobrevida, cohorte, thr)
  )

  print(curva_sobrevida)
  table_matrix <- table(stats_con_tiempo$group, stats_con_tiempo$evento)
  print(table_matrix)
}
