#Unificar fusion reports ---------------------------
library(OMICsdo)
RNAseqP("/media/4tb1/Daniela/Environ/Bulk-10-Colon", RunARRIBA = TRUE)

library(readxl)

Metadata_M <- as.data.frame(read_excel("/media/4tb1/Daniela/Environ/RedNeuronal/Metadata Environ -  Mama.xlsx"))
#Reemplazo los puntos por guiones:
Metadata_M$`ID ENVIRON` <- gsub("\\.", "-", Metadata_M$`ID ENVIRON`)
Metadata_C <- as.data.frame(read_excel("/media/4tb1/Daniela/Environ/RedNeuronal/Metadata Environ - Colon.xlsx"))

# Mama TCL -------------------------------------------------------------
patients_dir <- "/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/TCL_15/Muestras/Muestras_A/PrimeraTanda"
patient_dir <- "/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/TCL_15/Muestras/Muestras_A/PrimeraTanda/01/trimmed"
library(OMICsdo)
RNAseqP(patients_dir = patients_dir, RunARRIBA = TRUE)

#Mama CD GENOMICS ------------------------------------------------------
Todos_FusionReports_Mama <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/Todos-FusionReports-Mama.xlsx")
colnames(Todos_FusionReports_Mama)[1] <- "ID"
length(unique(Todos_FusionReports_Mama$ID))
Todos_FusionReports_Mama$Tecnologia <- "Illumina"
Todos_FusionReports_Mama$Secuenciador <- "CDGenomics"
Todos_FusionReports_Mama$ID[which(Todos_FusionReports_Mama$ID == "Br0307-2")] <- "Br0307"
Todos_FusionReports_Mama$ID[which(Todos_FusionReports_Mama$ID == "Br0310-1")] <- "Br0310"

Metadata_global <- read_excel("/media/16TBDisk/Daniela/Environ/Environ_Paula_metadata_global_2025.xlsx")
Metadata_global <- Metadata_global[which(Metadata_global$categoria_1 == "mama_61"),]
Metadata_global <- Metadata_global[which(Metadata_global$ID %in% Todos_FusionReports_Mama$ID),]
all(Todos_FusionReports_Mama$ID %in% Metadata_global$ID)
str(Metadata_global)
Metadata_global$tiempo <- as.numeric(Metadata_global$tiempo)

Todos_FusionReports_Mama <- merge(
  Todos_FusionReports_Mama,
  Metadata_global[, c("ID", "tiempo")],
  by = "ID",
  all.x = TRUE
)

Todos_FusionReports_Mama$tiempo

StatsFusions_Mama <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/StatsFusions-Mama.xlsx")
StatsFusions_Mama$Tecnologia <- "Illumina"
StatsFusions_Mama$Secuenciador <- "CDGenomics"
StatsFusions_Mama$ID[which(StatsFusions_Mama$ID == "Br0307-2")] <- "Br0307"
StatsFusions_Mama$ID[which(StatsFusions_Mama$ID == "Br0310-1")] <- "Br0310"
StatsFusions_Mama <- merge(StatsFusions_Mama,
                           Metadata_global[, c("ID", "tiempo", "evento", "G18G19")],
                           by = "ID",
                           all.x = TRUE)


#Analisis sobrevida:
any(is.na(StatsFusions_Mama$tiempo))
any(is.na(StatsFusions_Mama$evento))
StatsFusions_Mama <- StatsFusions_Mama[-which(is.na(StatsFusions_Mama$tiempo)),]

#solo teniendo en cuenta Her2, luminal/her2 y luminal B
unique(StatsFusions_Mama$Grupo)
StatsFusions_Mama <- StatsFusions_Mama[which(StatsFusions_Mama$Grupo %in% c("Luminal HER2","Luminal B", "HER2")),]

str(StatsFusions_Mama)
summary(StatsFusions_Mama$Fusiones_conf_H)
valor_corte <- 8
valor_corte <- 10.5

StatsFusions_Mama$group <- ifelse(StatsFusions_Mama$Fusiones_conf_H > valor_corte, "High", "Low")
StatsFusions_Mama$group <- factor(StatsFusions_Mama$group, levels = c("Low","High"))
#StatsFusions_Mama$V <- ifelse(StatsFusions_Mama$ == "Alive", 0, 1)

table(StatsFusions_Mama$group)

library("survival")
library("survminer")

#install.packages("/media/4tb2/Daniela/Environ/Fusiones/survminer_0.5.0.tar.gz", repos = NULL, type = "source")
#install.packages("/media/4tb2/Daniela/Environ/Fusiones/Matrix_1.6-0.tar.gz", repos = NULL, type = "source")
#install.packages(c("MatrixModels", "quantreg", "car", "rstatix", "ggpubr", "survminer"))

surv_object <- Surv(time = StatsFusions_Mama$tiempo, event = StatsFusions_Mama$evento)
fit1 <- survfit(surv_object ~ group, data = StatsFusions_Mama)
#ggsurvplot(fit1, data = StatsFusions_Mama, pval = TRUE, risk.table = TRUE)

ggsurvplot(fit1, data = StatsFusions_Mama, size = 1,  # change line size
           linetype = "strata", # change line type by groups
           palette = c("red","blue"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE # add table
)


table(StatsFusions_Mama$group, StatsFusions_Mama$MTT)
table(StatsFusions_Mama$MTT)

# Boxplot para comparar por grupos de MTT
ggplot(StatsFusions_Mama, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_boxplot() +
  labs(
    title = "TFB por MTT",
    x = "Metástasis",
    y = "Fusiones con HIGH confidence"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) # Centrar el título

wilcoxon_test_result <- wilcox.test(Fusiones_conf_H ~ MTT, data = StatsFusions_Mama)
wilcoxon_test_result$p.value

# Boxplot para comparar por Subtipo
ggplot(StatsFusions_Mama, aes(x = Grupo, y = Fusiones_conf_H, fill = Grupo)) +
  geom_boxplot() +
  labs(
    title = "TFB por Subtipo",
    x = "Subtipo",
    y = "Fusiones con HIGH confidence"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) # Centrar el título


kruskal_test_result <- kruskal.test(Fusiones_conf_H ~ Grupo, data = StatsFusions_Mama)
kruskal_test_result$p.value



#--------------------------------------------------------------

patients_dir <- "/media/16TBDisk/Daniela/Environ/TodoMama"
ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
#id = ids[2]
length(ids)

#Colon
patients_dir <- "/media/16TBDisk/Daniela/Environ/TodoColon"
ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
ids2 <- list.dirs(path = "/media/4tb1/Daniela/Environ/Bulk-10-Colon", full.names = TRUE, recursive = FALSE)
length(ids2)
ids <- c(ids, ids2)
length(ids)


#SRA
patients_dir <- "/media/4tb1/Daniela/Environ/Fusiones/MuestrasSRA"
patients_dir <- "/media/4tb2/Daniela/MuestrasChile2/Muestras"
ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
length(ids)
ids <- ids[-1]

MetadataSRA <- read.table("/media/4tb1/Daniela/Environ/Fusiones/SraRunTable.txt", header = TRUE, sep = ",")
colnames(MetadataSRA)[which(colnames(MetadataSRA) == "metastasis")] <- "MTT"
colnames(MetadataSRA)[which(colnames(MetadataSRA) == "Run")] <- "ID"

Todos_FusionReport <- data.frame()
Stats_Fusions <- data.frame()
k=1
#id <- ids[2]
for (id in ids) {
  i <- basename(id)
  print(i)
  fusions_file <- sprintf("%s/trimmed/%s_FusionReport.xlsx", id, i)
  FusionReport <- read_excel(fusions_file)

  if(!(nrow(FusionReport) == 0)) {
    if (grepl("Co", i)) {
      Metadata <- Metadata_C
      Grupo <- Metadata[which(Metadata$ID == i), "Grupo"]
    } else if (grepl("Br", i)) {
      Metadata <- Metadata_M
      Grupo <- Metadata[which(Metadata$ID == i), "Subtipo"]
    } else {
      Metadata <- MetadataSRA
      Grupo <- Metadata[which(Metadata$ID == i), "tissue"]
    }

    FusionReport$MTT <- Metadata[which(Metadata$ID == i), "MTT"]
    FusionReport$Grupo <- Grupo
    FusionReport$ID <-i

    FusionReport_ID <- cbind(FusionReport$ID, FusionReport[,1:(ncol(FusionReport)-1)])
    colnames(FusionReport)[1] <- "ID"

    if(nrow(Todos_FusionReport) == 0){
      Todos_FusionReport <- FusionReport_ID
    } else {
      Todos_FusionReport <- rbind(Todos_FusionReport, FusionReport_ID)
    }

  } else {
    print(sprintf("Se encontraron 0 fusiones en %s", i))
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


  if (grepl("Co", i)) {
    Metadata <- Metadata_C
    Grupo <- Metadata[which(Metadata$ID == i), "Grupo"]
  } else if (grepl("Br", i)){
    Metadata <- Metadata_M
    Grupo <- Metadata[which(Metadata$ID == i), "Subtipo"]
  } else {
    Metadata <- MetadataSRA
    Grupo <- Metadata[which(Metadata$ID == i), "tissue"]
  }

  colnames(Metadata)[1] <- "ID"
  Stats_Fusions[k, "MTT"] <- Metadata[which(Metadata$ID == i), "MTT"]
  Stats_Fusions[k, "Grupo"] <- Grupo

  k = k+1

}

write.xlsx(Todos_FusionReport, file = "/media/16TBDisk/Daniela/Environ/TodoColon/Todos-FusionReports-Colon.xlsx")
write.xlsx(Stats_Fusions, file = "/media/16TBDisk/Daniela/Environ/TodoColon/StatsFusions-Colon.xlsx")

write.xlsx(Todos_FusionReport, file = "/media/16TBDisk/Daniela/Environ/TodoMama/Todos-FusionReports-Mama.xlsx")
write.xlsx(Stats_Fusions, file = "/media/16TBDisk/Daniela/Environ/TodoMama/StatsFusions-Mama.xlsx")

write.xlsx(Todos_FusionReport, file = "/media/4tb1/Daniela/Environ/Fusiones/FusionReports-52SRA.xlsx")
write.xlsx(Stats_Fusions, file = "/media/4tb1/Daniela/Environ/Fusiones/StatsFusions-52SRA.xlsx")

library(googledrive)
#drive_upload("/media/16TBDisk/Daniela/Environ/TodoColon/StatsFusions-Colon.xlsx", name = "StatsFusions-Colon.xlsx")


summary <- data.frame(do.call(cbind, lapply(Stats_Fusions, summary)))

length(unique(Todos_FusionReport$`FusionReport$ID`))
length(unique(Todos_FusionReport$`FusionReport$ID`[which(Todos_FusionReport$confidence == "high")]))
muestras_con_high <- unique(Todos_FusionReport$`FusionReport$ID`[which(Todos_FusionReport$confidence == "high")])
muestras_sin_high <- unique(Todos_FusionReport$`FusionReport$ID`[-which(Todos_FusionReport$`FusionReport$ID`%in% muestras_con_high)])
info_muestras_sin_high <- Todos_FusionReport[which(Todos_FusionReport$`FusionReport$ID` %in% muestras_sin_high), c(1,2,3,6,7, 16, 38, 37)]

write.xlsx(info_muestras_sin_high, file = "/media/16TBDisk/Daniela/Environ/TodoMama/Info_Muestras_SinFusionesHIGH.xlsx")
write.xlsx(info_muestras_sin_high, file = "/media/16TBDisk/Daniela/Environ/TodoColon/Info_Muestras_SinFusionesHIGH.xlsx")


# Crear un boxplot para cada grupo, comparando MTT+ y MTT-
library(ggplot2)

# FUSIONES HIGH -------------------------------------
#Grupos originales:
library(ggplot2)
ggplot(Stats_Fusions, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +  # Gráfico de violín con un poco de transparencia
  #geom_boxplot(width = 0.1, position = position_dodge(0.75), outlier.shape = NA) +  # Boxplot sin los outliers para no sobrecargar el gráfico
  geom_boxplot()+
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Puntos para cada muestra con algo de dispersión y transparencia
  #facet_wrap(~ Grupo, scales = "free_y") +
  facet_wrap(~ Grupo) +
  labs(title = "Distribución de la Cantidad de Fusiones HIGH por Grupo y MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal()

ggplot(Stats_Fusions, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +  # Gráfico de violín con un poco de transparencia
  #geom_boxplot(width = 0.1, position = position_dodge(0.75), outlier.shape = NA) +  # Boxplot sin los outliers para no sobrecargar el gráfico
  geom_boxplot()+
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Puntos para cada muestra con algo de dispersión y transparencia
  #facet_wrap(~ Grupo, scales = "free_y") +
  labs(title = "Distribución de la Cantidad de Fusiones HIGH -  MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal()

wilcox.test(Fusiones_conf_H ~ MTT, data = Stats_Fusions)

unique_groups <- unique(Stats_Fusions$Grupo)
for (group in unique_groups) {
  group_data <- subset(Stats_Fusions, Grupo == group)
  test_result <- wilcox.test(Fusiones_conf_H ~ MTT, data = group_data)
  print(group)
  print(test_result)
}

# FUSIONES TOTALES  -------------------------------------
#Grupos originales:
ggplot(Stats_Fusions, aes(x = MTT, y = Cantidad_Fusiones, fill = MTT)) +
  geom_violin(alpha = 0.5) +  # Gráfico de violín con un poco de transparencia
  #geom_boxplot(width = 0.1, position = position_dodge(0.75), outlier.shape = NA) +  # Boxplot sin los outliers para no sobrecargar el gráfico
  geom_boxplot()+
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Puntos para cada muestra con algo de dispersión y transparencia
  facet_wrap(~ Grupo) +
  labs(title = "Distribución de la Cantidad de Fusiones TOTALES por Grupo y MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal()

ggplot(Stats_Fusions, aes(x = MTT, y = Cantidad_Fusiones, fill = MTT)) +
  geom_violin(alpha = 0.5) +  # Gráfico de violín con un poco de transparencia
  #geom_boxplot(width = 0.1, position = position_dodge(0.75), outlier.shape = NA) +  # Boxplot sin los outliers para no sobrecargar el gráfico
  geom_boxplot()+
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Puntos para cada muestra con algo de dispersión y transparencia
  #facet_wrap(~ Grupo, scales = "free_y") +
  labs(title = "Distribución de la Cantidad de Fusiones TOTALES por Grupo y MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal()
wilcox.test(Cantidad_Fusiones ~ MTT, data = Stats_Fusions)

unique_groups <- unique(Stats_Fusions$Grupo)
for (group in unique_groups) {
  group_data <- subset(Stats_Fusions, Grupo == group)
  test_result <- wilcox.test(Cantidad_Fusiones ~ MTT, data = group_data)
  print(group)
  print(test_result)
}

#Mama - 4 grupos -----------------------------------
Stats_Fusions$Grupo2 <- ifelse(Stats_Fusions$Grupo == "Luminal HER2", "HER2", Stats_Fusions$Grupo)
unique(Stats_Fusions$Grupo2 )
Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 == "HER2")] <- "Luminal HER2/HER2"

ggplot(Stats_Fusions, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +  # Gráfico de violín con un poco de transparencia
  #geom_boxplot(width = 0.1, position = position_dodge(0.75), outlier.shape = NA) +  # Boxplot sin los outliers para no sobrecargar el gráfico
  geom_boxplot()+
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Puntos para cada muestra con algo de dispersión y transparencia
  facet_wrap(~ Grupo2, scales = "free_y") +
  labs(title = "Distribución de la Cantidad de Fusiones HIGH por Grupo y MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal()
unique_groups <- unique(Stats_Fusions$Grupo2)
for (group in unique_groups) {
  group_data <- subset(Stats_Fusions, Grupo2 == group)
  test_result <- wilcox.test(Fusiones_conf_H ~ MTT, data = group_data)
  print(group)
  print(test_result)
}

#################################

boxplot(Stats_Fusions$Cantidad_Fusiones ~ Stats_Fusions$MTT)
wilcox.test(Stats_Fusions$Cantidad_Fusiones ~ Stats_Fusions$MTT)

boxplot(Stats_Fusions$Cantidad_Fusiones ~ Stats_Fusions$Grupo)
kruskal.test(Stats_Fusions$Cantidad_Fusiones ~ Stats_Fusions$Grupo)
boxplot(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo %in% c("HER2", "Luminal A"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("HER2", "Luminal A"))])
wilcox.test(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo %in% c("HER2", "Luminal A"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("HER2", "Luminal A"))])

boxplot(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$Grupo)
kruskal.test(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$Grupo)


boxplot(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$Grupo2)
kruskal.test(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$Grupo2)

boxplot(Stats_Fusions$Cantidad_Fusiones ~ Stats_Fusions$Grupo2)
kruskal.test(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$Grupo2)

boxplot(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))])
wilcox.test(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))])

boxplot(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))])
wilcox.test(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))])

boxplot(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))])
wilcox.test(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))])

#Para high:
boxplot(Stats_Fusions$Fusiones_conf_H[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))])
wilcox.test(Stats_Fusions$Fusiones_conf_H[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Triple Negativo"))])

boxplot(Stats_Fusions$Fusiones_conf_H[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))])
wilcox.test(Stats_Fusions$Fusiones_conf_H[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal A/B", "Luminal HER2/HER2"))])

boxplot(Stats_Fusions$Fusiones_conf_H[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))])
wilcox.test(Stats_Fusions$Fusiones_conf_H[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))] ~ Stats_Fusions$Grupo2[which(Stats_Fusions$Grupo2 %in% c("Luminal HER2/HER2", "Triple Negativo"))])

boxplot(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$MTT)
wilcox.test(Stats_Fusions$Fusiones_conf_H ~ Stats_Fusions$MTT)

boxplot(Stats_Fusions$Fusiones_conf_M ~ Stats_Fusions$MTT)
wilcox.test(Stats_Fusions$Fusiones_conf_M ~ Stats_Fusions$MTT)
boxplot(Stats_Fusions$Fusiones_conf_L ~ Stats_Fusions$MTT)
wilcox.test(Stats_Fusions$Fusiones_conf_L ~ Stats_Fusions$MTT)


#Para colon
#total de fusiones --
boxplot(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))])
wilcox.test(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))])
boxplot(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo III"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo III"))])
wilcox.test(Stats_Fusions$Cantidad_Fusiones[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo III"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo III"))])

boxplot(Stats_Fusions$Fusiones_conf_M[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))])
wilcox.test(Stats_Fusions$Fusiones_conf_M[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo I", "Grupo IV"))])
boxplot(Stats_Fusions$Fusiones_conf_M[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo IV"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo IV"))])
wilcox.test(Stats_Fusions$Fusiones_conf_M[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo IV"))] ~ Stats_Fusions$Grupo[which(Stats_Fusions$Grupo %in% c("Grupo II", "Grupo IV"))])



# Frecuencia de fusiones  ------------------------------------------------------------
#Todas las fusiones
algunas <- Todos_FusionReport[, c(1,2,3,6,7, 38, 37)]
colnames(algunas)
library(dplyr)

# Seleccionar las columnas de interés y contar las repeticiones de filas
Freq_fusiones <- algunas %>%
  select(gene1, gene2, breakpoint1, breakpoint2) %>%
  group_by(across(everything())) %>%
  summarise(count = n(), .groups = 'drop')


#Fusiones con high confidence
algunas <- Todos_FusionReport[which(Todos_FusionReport$confidence == "high"), c(1,2,3,6,7, 16, 38, 37)]
colnames(algunas)
class(algunas)
install.packages("dplyr")
library(dplyr)

# Seleccionar las columnas de interés y contar las repeticiones de filas
Freq_fusiones_high <- algunas %>%
  select(gene1, gene2, breakpoint1, breakpoint2) %>%
  group_by(across(everything())) %>%
  summarise(count = n(), .groups = 'drop')

colnames(Freq_fusiones_high)[ncol(Freq_fusiones_high)] <- "Freq_FusionesHigh"

library(openxlsx)
write.xlsx(Freq_fusiones_high, file = "/media/16TBDisk/Daniela/Environ/TodoMama/Frecuencia_FusionesHigh-Mama.xlsx")
write.xlsx(Freq_fusiones_high, file = "/media/16TBDisk/Daniela/Environ/TodoColon/Frecuencia_FusionesHigh-Colon.xlsx")
write.xlsx(Freq_fusiones_high, file = "/media/4tb1/Daniela/Environ/Fusiones/Frecuencia_FusionesHigh-SRA.xlsx")


#Frecuencia de genes ------------------------------------------------------------------
algunas <- Todos_FusionReport[which(Todos_FusionReport$confidence == "high"), c(1,2,3,6,7, 16, 38, 37)]
genes <- c(algunas$gene1, algunas$gene2)
length(genes)
split_genes <- unlist(strsplit(genes, split = ","))
length(split_genes)
# eliminar la información entre paréntesis
gene_names <- gsub("\\(.*?\\)", "", split_genes)
# Quitar espacios adicionales si existen
gene_names <- trimws(gene_names)

length(unique(gene_names))
nrow(Todos_FusionReport)

#Cantidad de fusiones en las que participa cada gen:
Frec_g <- as.data.frame(table(gene_names))
Frec_g$PorcentajeFusiones<- round((Frec_g$Freq/nrow(algunas))*100, 2)

#Cantidad de muestras en las que aparece cada gen: -----
Frec_g$Cantidad_Muestras <- 0
Frec_g$Cantidad_Fusiones <- 0
Frec_g$MTTneg <-0
Frec_g$MTTpos <-0

Frec_g$Normal <-0
Frec_g$Tumor <-0

Frec_g$LuminalA <-0
Frec_g$LuminalB <-0
Frec_g$LuminalHER2_HER2 <-0
Frec_g$TripleNegativo <-0

Frec_g$GrupoI <- 0
Frec_g$GrupoII <- 0
Frec_g$GrupoIII <- 0
Frec_g$GrupoIV <- 0

i=1
for (i in 1:nrow(Frec_g)){
  gen <- Frec_g[i,1]
  print(gen)

  #Tengo que convertir las celdas gene1 y gene2 en vectores con los nombres de los genes separados:

  k=1
  filas_si_dfgen <- c()
  for(k in 1:nrow(algunas)) {
    genes_split <- unlist(strsplit(algunas$gene1[k], split = ","))
    # eliminar la información entre paréntesis
    genes_split <- gsub("\\(.*?\\)", "", genes_split)
    # Quitar espacios adicionales si existen
    genes_split <- trimws(genes_split)

    if (any(gen == genes_split)) {
      filas_si_dfgen <- c(filas_si_dfgen, k)
    }
  }

  for(k in 1:nrow(algunas)) {
    genes_split <- unlist(strsplit(algunas$gene2[k], split = ","))
    # eliminar la información entre paréntesis
    genes_split <- gsub("\\(.*?\\)", "", genes_split)
    # Quitar espacios adicionales si existen
    genes_split <- trimws(genes_split)

    if (any(gen == genes_split)) {
      filas_si_dfgen <- c(filas_si_dfgen, k)
    }
  }

  df_gen <- algunas[filas_si_dfgen,]

  cant_muestras_unicas <- length(unique(df_gen$`FusionReport$ID`))
  muestras_unicas <- unique(df_gen$`FusionReport$ID`)

  MTT_muestras_unicas <- c()
  Grupo_muestras_unicas <- c()
  #muestra_unica <- muestras_unicas[1]
  for (muestra_unica in muestras_unicas){
    print(muestra_unica)
    MTT <- unique(df_gen$MTT[which(df_gen$`FusionReport$ID` == muestra_unica)])
    print(MTT)
    MTT_muestras_unicas <- c(MTT_muestras_unicas, MTT)
    Grupo <- unique(df_gen$Grupo[which(df_gen$`FusionReport$ID` == muestra_unica)])
    Grupo_muestras_unicas <- c(Grupo_muestras_unicas, Grupo)
  }

  Frec_g$Cantidad_Muestras[i] <- cant_muestras_unicas
  Frec_g$PorcentajeMuestras[i] <- (cant_muestras_unicas/length(unique(algunas$`FusionReport$ID`)))*100
  Frec_g$Cantidad_Fusiones[i] <- nrow(df_gen)

  contadorMTTneg <- 0
  contadorMTTpos <- 0
  for (cada_mtt in MTT_muestras_unicas) {
    if(cada_mtt == "MET-" | cada_mtt == "no") {
      contadorMTTneg <- contadorMTTneg + 1
    } else if(cada_mtt == "MET+" | cada_mtt == "yes") {
      contadorMTTpos <- contadorMTTpos + 1
    }
  }
  Frec_g$MTTneg[i] <- contadorMTTneg
  Frec_g$MTTpos[i] <- contadorMTTpos

  if (grepl("Br", algunas$`FusionReport$ID`[i])) {
    contadorLA <- 0
    contadorLB <- 0
    contadorHER2 <- 0
    contadorTN <- 0
    for(cada_grupo in Grupo_muestras_unicas) {
      if(cada_grupo == "Luminal A") {
        contadorLA <- contadorLA + 1
      } else if(cada_grupo == "Luminal B") {
        contadorLB <- contadorLB + 1
      } else if(cada_grupo == "Luminal HER2" | cada_grupo == "HER2") {
        contadorHER2 <- contadorHER2 + 1
      } else if (cada_grupo == "Triple Negativo") {
        contadorTN <- contadorTN + 1
      }
    }
    Frec_g$LuminalA[i] <- contadorLA
    Frec_g$LuminalB[i] <- contadorLB
    Frec_g$LuminalHER2_HER2[i] <- contadorHER2
    Frec_g$TripleNegativo[i] <- contadorTN

  } else if (grepl("Co", algunas$`FusionReport$ID`[i])){
    contadorI <- 0
    contadorII <- 0
    contadorIII <- 0
    contadorIV <- 0
    for(cada_grupo in Grupo_muestras_unicas) {
      if(cada_grupo == "Grupo I") {
        contadorI <- contadorI + 1
      } else if(cada_grupo == "Grupo II") {
        contadorII <- contadorII + 1
      } else if(cada_grupo == "Grupo III") {
        contadorIII <- contadorIII + 1
      } else if (cada_grupo == "Grupo IV") {
        contadorIV <- contadorIV + 1
      }
    }
    Frec_g$GrupoI[i] <- contadorI
    Frec_g$GrupoII[i] <- contadorII
    Frec_g$GrupoIII[i] <- contadorIII
    Frec_g$GrupoIV[i] <- contadorIV

  } else {
    contadorN <- 0
    contadorT <- 0

    for(cada_grupo in Grupo_muestras_unicas) {
      if(cada_grupo == "breast tumor") {
        contadorT <- contadorT + 1
      } else if(cada_grupo == "normal breast tissue") {
        contadorN <- contadorN + 1
      }
    }
    Frec_g$Normal[i] <- contadorN
    Frec_g$Tumor[i] <- contadorT
  }

}


# EL PROBLEMA ESTÁ EN DETECTAR LOS NOMBRES DE LOS GENES EN EL VECTOR SEPARADO POR COMAS!!!!!!!!!!!1
all(Frec_g$Freq == Frec_g$Cantidad_Fusiones)
f_raro <- Frec_g[which(!(Frec_g$Freq == Frec_g$Cantidad_Fusiones)),]

write.xlsx(Frec_g, file = "/media/16TBDisk/Daniela/Environ/TodoMama/GenesFusionesHIGH-Mama.xlsx")
write.xlsx(Frec_g, file = "/media/16TBDisk/Daniela/Environ/TodoColon/GenesFusionesHIGH-Colon.xlsx")
write.xlsx(Frec_g, file = "/media/4tb1/Daniela/Environ/Fusiones/GenesFusionesHIGH-SRA.xlsx")

library(readxl)
GenesFusionesH <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/GenesFusionesHIGH-Mama.xlsx")
GenesFusionesHIGH_SRA <- read_excel("/media/4tb1/Daniela/Environ/Fusiones/GenesFusionesHIGH-SRA.xlsx")

genes_comun <- intersect(GenesFusionesH$gene_names, Frec_g$gene_names)
Frec_comunEnviron <- GenesFusionesH[which(GenesFusionesH$gene_names %in% genes_comun), c(1,4)]
Frec_comunSRA <- Frec_g[which(Frec_g$gene_names %in% genes_comun), c(1,4)]

ggplot(GenesFusionesH, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +  # Gráfico de violín con un poco de transparencia
  #geom_boxplot(width = 0.1, position = position_dodge(0.75), outlier.shape = NA) +  # Boxplot sin los outliers para no sobrecargar el gráfico
  geom_boxplot()+
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +  # Puntos para cada muestra con algo de dispersión y transparencia
  facet_wrap(~ Grupo2, scales = "free_y") +
  labs(title = "Distribución de la Cantidad de Fusiones HIGH por Grupo y MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal()


library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")

# Chart
venn.plot <- venn.diagram(
  x = list(GenesFusionesH$gene_names, Frec_g$gene_names),
  category.names = c("Environ" , "SRA"),
  filename = NULL,
  output=TRUE,
  # Output features
  imagetype="png" ,
  height = 20 ,
  width = 20 ,
  resolution = 300,
  compression = "lzw",

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],

  # Numbers
  cex = 2
  #fontface = "bold",
  #fontfamily = "sans"

  # Set names
)

grid.draw(venn.plot)



library(reshape2)
solo_MTT <- Frec_g[, c(1,6,7)]
df_long <- melt(solo_MTT, id.vars = "gene_names", variable.name = "Condition", value.name = "Count")

solo_Subtipo <- Frec_g[, c(1,8,9,10,11)]
df_long <- melt(solo_Subtipo, id.vars = "gene_names", variable.name = "Condition", value.name = "Count")

# Crear el boxplot
ggplot(df_long, aes(x = Condition, y = Count, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Boxplot de Cantidad de Fusiones en Muestras MTT Pos y Neg",
       x = "Condición", y = "Cantidad de Fusiones") +
  theme_minimal()

ggplot(df_long, aes(x = Condition, y = Count, fill = Condition)) +
  geom_boxplot() +
  labs(title = "Boxplot de Cantidad de Fusiones en Muestras según Subtipos",
       x = "Condición", y = "Cantidad de Fusiones") +
  theme_minimal()

#########################################################
# Fusiones en dataset publico
#https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA762469&o=acc_s%3Aa

#Descarga metadata: -----------------------
#https://www.ncbi.nlm.nih.gov/Traces/study/?page=2&acc=PRJNA762469&o=metastasis_sam_s%3Aa%3Bacc_s%3Aa
sra_metadata <- read.delim("/media/4tb1/Daniela/Environ/Fusiones/SraRunTable.txt", header = TRUE, sep = ",")


# Descarga desde TXT de SRA: ---------------------
options(timeout = 100000)
sra_urls <- read.delim("/media/4tb1/Daniela/Environ/Fusiones/sra_explorer_fastq_urls - NUEVA.txt", header = FALSE)
sra_urls <- as.data.frame(sra_urls[-c(1,2,4,5,6,7,8,9),1])

i=80
for (i in 80:nrow(sra_urls)) {
#for (i in 1:10) {
  fastq <- sra_urls[i,]
  id <- strsplit(basename(fastq), split="_")[[1]][1]
  R <- ifelse(grepl(1, strsplit(basename(fastq), split="_")[[1]][2]), "R1", "R2")
  dir.create(sprintf("/media/4tb2/Daniela/MuestrasSRA-Nuevo/%s", id))

  if(!(file.exists(sprintf("/media/4tb2/Daniela/MuestrasSRA-Nuevo/%s/%s_%s.fastq.gz", id, id, R)))) {
    download.file(fastq, destfile = sprintf("/media/4tb2/Daniela/MuestrasSRA-Nuevo/%s/%s_%s.fastq.gz", id, id, R))
    setwd("/media/4tb2/Daniela/MuestrasSRA-Nuevo")
    #system("wget", fastq, intern=FALSE, ignore.stdout = TRUE, ignore.stderr= TRUE, wait= TRUE)
    #system(sprintf("wget -O %s %s", destfile, fastq))
  } else {
    print(sprintf("La muestra %s %s ya fue descargada", id, R))
  }
  i <- i+1
}

library(OMICsdo)
patient_dir <- "/media/4tb1/Daniela/Environ/Fusiones/MuestrasSRA-Nuevo/SRR15852394/trimmed"
RNAseqP("/media/4tb1/Daniela/Environ/Fusiones/MuestrasSRA", fastQC_after_trim = FALSE, plot_FastQC_trim = FALSE, RunARRIBA = TRUE)
RNAseqP("/media/4tb1/Daniela/Environ/Fusiones/Problematicos", fastQC_after_trim = FALSE, plot_FastQC_trim = FALSE, RunARRIBA = TRUE)

#SRA
patients_dir <- "/media/4tb1/Daniela/Environ/Fusiones/MuestrasSRA"
ids <-  list.dirs(path = patients_dir, full.names = TRUE, recursive = FALSE)
length(ids)

Todos_FusionReport <- data.frame()
Stats_Fusions <- data.frame()
k=1
id <- ids[2]
for (id in ids) {
  i <- basename(id)
  print(i)
  fusions_file <- sprintf("%s/trimmed/%s_FusionReport.xlsx", id, i)
  FusionReport <- read_excel(fusions_file)

  FusionReport$MTT <- sra_metadata[which(sra_metadata$Run == i), "metastasis"]
  FusionReport$Tejido <- sra_metadata[which(sra_metadata$Run == i), "tissue"]
  FusionReport$ID <-i

  FusionReport_ID <- cbind(FusionReport$ID, FusionReport[,1:(ncol(FusionReport)-1)])
  colnames(FusionReport)[1] <- "ID"

  if(nrow(Todos_FusionReport) == 0){
    Todos_FusionReport <- FusionReport_ID
  } else {
    Todos_FusionReport <- rbind(Todos_FusionReport, FusionReport_ID)
  }

  #Estadisticas de los reportes -----------------------------------

  Stats_Fusions[k, "ID"] <- i
  Stats_Fusions[k, "Cantidad_Fusiones"] <- nrow(FusionReport)

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

  k = k+1

}

GenesFusionesHIGH_SRA <- read_excel("/media/4tb1/Daniela/Environ/Fusiones/GenesFusionesHIGH-SRA.xlsx")
GenesFusionesHIGH_SRA <- read_excel("/media/4tb1/Daniela/Environ/Fusiones/GenesFusionesHIGH-SRA.xlsx")

GenesFusionesH <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/GenesFusionesHIGH-Mama.xlsx")
Todos_FusionReport <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/Todos-FusionReports-Mama.xlsx")
Stats_Fusions <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/StatsFusions-Mama.xlsx")
