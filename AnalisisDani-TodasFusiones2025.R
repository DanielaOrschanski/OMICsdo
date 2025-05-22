library(OMICsdo)
RNAseqP(patients_dir = "/media/4tb1/Paula/13Br_Validacion/Muestras", RunARRIBA = TRUE)

list_p <- list.dirs("/media/4tb1/Paula/13Br_Validacion/Muestras", recursive = FALSE, full.names = TRUE)
ids <- list.dirs("/media/4tb1/Paula/13Br_Validacion/Muestras", recursive = FALSE, full.names = FALSE)
#p <- list_p[2]
for(p in list_p) {
  print(p)
  ptrim <- paste0(p, "/trimmed", sep ="")
  runARRIBA(ptrim)
}

Environ_Paula_metadata_global_2025 <- read_excel("/media/16TBDisk/Daniela/Environ/Environ_Paula_metadata_global_2025.xlsx")

# CDGenomics --------------------------------------
metadata_CDG <- Environ_Paula_metadata_global_2025[which(Environ_Paula_metadata_global_2025$categoria_1 == "val_20"),]
ids <- list.dirs("/media/4tb1/Paula/13Br_Validacion/Muestras", recursive = FALSE, full.names = FALSE)
all( ids %in% metadata_CDG$ID)
out_CGD <- fusionStats(patients_dir = "/media/4tb1/Paula/13Br_Validacion/Muestras",
                       Metadata = metadata_CDG,
                       group ="subtipo",
                       cohorte = "CDGenomics",
                       sobrevida = FALSE)
fusions_CGD <- out_CGD[[1]]
stats_CGD <- out_CGD[[2]]
genes_CGD <- out_CGD[[3]]

write.xlsx(TFB_china, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_TFB.xlsx")
write.xlsx(metadata_china_completa, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_Metadata.xlsx")
write.xlsx(genes_china, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_Genes.xlsx")
write.xlsx(stats_china, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_FusionStats.xlsx")


# ------------------------

metadata_mama61 <- Environ_Paula_metadata_global_2025[which(Environ_Paula_metadata_global_2025$categoria_1 == "mama_61"),]
colnames(metadata_mama61)
metadata_mama61$ID[which(metadata_mama61$ID == "Br0310")] <- "Br0310-1"
metadata_mama61$ID[which(metadata_mama61$ID == "Br0307")] <- "Br0307-2"
library(ggpubr)
out_mama61 <- fusionStats(patients_dir = "/media/16TBDisk/Daniela/Environ/TodoMama",
                          Metadata = metadata_mama61,
                          group = "subtipo",
                          cohorte = "Mama61",
                          sobrevida = TRUE)
stats_mama61 <- out_mama61[[2]]
gen_mama61 <- out_mama61[[3]]
fusions_mama61 <- out_mama61[[1]]

length(which(grepl("kinase", fusions_mama61$retained_protein_domains))) / nrow(fusions_mama61) *100
unique(fusions_mama61$confidence)
length(which(grepl("kinase", fusions_mama61$retained_protein_domains[fusions_mama61$confidence == "high"]))) / nrow(fusions_mama61[fusions_mama61$confidence == "high",]) *100

write.csv(TFB_mama61, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Environ_Mama61_TFB.csv")
write.xlsx(TFB_mama61, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Environ_Mama61_TFB.xlsx")
write.xlsx(metadata_mama61, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Environ_Mama61_Metadata.xlsx")
write.xlsx(gen_mama61, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Environ_Mama61_Genes.xlsx")
write.xlsx(stats_mama61, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Environ_Mama61_FusionStats.xlsx")

TFB_mama61 <- fusions_mama61[fusions_mama61$confidence == "high",]

library(dplyr)
colnames(stats_mama61)

conteo_mtt_por_grupo <- stats_mama61 %>%
  group_by(Grupo, MTT) %>%
  summarise(n = n(), .groups = "drop")

print(conteo_mtt_por_grupo)

summary_stats <- stats_mama61 %>%
  #group_by(MTT) %>%
  group_by(Grupo, MTT) %>%
  summarise(
    n_pacientes = n(),
    fusiones_media = mean(Fusiones_conf_H, na.rm = TRUE),
    fusiones_mediana = median(Fusiones_conf_H, na.rm = TRUE),
    #fusiones_sd = sd(Fusiones_conf_H, na.rm = TRUE),
    fusiones_min = min(Fusiones_conf_H, na.rm = TRUE),
    fusiones_q1 = quantile(Fusiones_conf_H, 0.25, na.rm = TRUE),
    fusiones_q3 = quantile(Fusiones_conf_H, 0.75, na.rm = TRUE),
    fusiones_max = max(Fusiones_conf_H, na.rm = TRUE)
  )

print(summary_stats)

stats_mama61$Predicho <- ifelse(stats_mama61$Fusiones_conf_H >= 8, "MET+", "MET-")
library(caret)
conf_matrix <- confusionMatrix(
  table(stats_mama61$Predicho, stats_mama61$MTT)
)
print(conf_matrix)
ggplot(stats_mama61, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "Cantidad de Fusiones HIGH por MTT - MAMA61",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal() +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 5, 10, 20, 30))

#-----------------------------------------------------------------------------------------------------
RNAseqP(patients_dir = "/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/Muestras",
        RunARRIBA = TRUE)

metadata_tcl <- Environ_Paula_metadata_global_2025[which(Environ_Paula_metadata_global_2025$categoria_1 == "tcl"),]
ids <- list.dirs("/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/Muestras", recursive = FALSE, full.names = FALSE)
metadata_tcl$ID <- gsub("TCL_", "", metadata_tcl$ID)
all( ids %in% metadata_tcl$ID)
out_tcl <- fusionStats(patients_dir = "/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/Muestras",
                       Metadata = metadata_tcl,
                       group = NA,
                       cohorte = "TCL",
                       sobrevida = FALSE)
stats_tcl <- out_tcl[[2]]
#stats_tcl <- read_excel("/media/4tb1/Paula/28Br_Validacion_TCL_Nuevas/Muestras/StatsFusions.xlsx")
stats_tcl$Fusiones_conf_H_log <- log10(stats_tcl$Fusiones_conf_H + 1)

ggplot(stats_tcl, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "Cantidad de Fusiones HIGH por MTT - TCL",
       x = "MTT",
       y = "Cantidad de Fusiones (pseudolog)") +
  theme_minimal() +
  #scale_y_log10()
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 5, 10, 20, 30))
  #scale_y_continuous(limits = c(0, NA))


#--------------------
library(readr)
BRCA_TCGA <- read_csv("/media/4tb1/Paula/BRCA_TCGA.csv")
metadata_tcga <- readRDS("/media/4tb1/Paula/BRCA_rna.rds")
metadata_tcga <- metadata_tcga$targets
metadata_tcga$barcode <-  gsub("\\.", "-", metadata_tcga$barcode)
colnames(BRCA_TCGA)[2] <- "ID"
colnames(BRCA_TCGA)[3] <- "Fusiones_conf_H"
colnames(BRCA_TCGA)[which(colnames(BRCA_TCGA) == "ajcc_pathologic_m")] <- "MTT"
BRCA_TCGA <- BRCA_TCGA[which(BRCA_TCGA$MTT %in% c("M0", "M1")),]
table(BRCA_TCGA$MTT)
BRCA_TCGA$MTT <- ifelse(BRCA_TCGA$MTT == "M0", "MET-", "MET+")
all(BRCA_TCGA$ID %in% metadata_tcga$barcode)
colnames(metadata_tcga)[1] <- "ID"
BRCA_TCGA <- merge(BRCA_TCGA, metadata_tcga[,c("ID", "pam50")], by= "ID")

stats_TCGA <- BRCA_TCGA[, c("ID", "Fusiones_conf_H", "MTT", "pam50")]
str(stats_TCGA)
colnames(stats_TCGA)[4] <- "Grupo"
unique(stats_TCGA$Grupo)
stats_TCGA <- stats_TCGA[which(stats_TCGA$Grupo %in% c("Her2", "LumB")),]
stats_TCGA <- stats_TCGA[stats_TCGA$MTT == "MET+" |
                           (stats_TCGA$MTT == "MET-" & stats_TCGA$Grupo %in% c("Her2", "LumB")), ]

table(stats_TCGA$MTT)
table(stats_TCGA$Grupo)
library(dplyr)

write.xlsx(metadata_tcga, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_TCGA_Metadata.xlsx")
write.xlsx(stats_TCGA, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_TCGA_FusionStats.xlsx")


conteo_mtt_por_grupo <- stats_TCGA %>%
  group_by(Grupo, MTT) %>%
  summarise(n = n(), .groups = "drop")

print(conteo_mtt_por_grupo)

ggplot(stats_TCGA, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "Cantidad de Fusiones HIGH por MTT - TCGA",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))+
  facet_grid( ~ Grupo)

stats_52SRA <- read_excel("/media/4tb2/Daniela/Environ/FusionesActualizado/StatsFusions-52SRA.xlsx")


# nottingham ---
metadata_nott <- read_excel("/media/4tb1/Paula/33Br_Nottingham/metadata_curada.xlsx")
colnames(metadata_nott)[3] <- "ID"
colnames(metadata_nott)[11] <- "MTT"
setwd(patients_dir)

stats_nottingham <- as.data.frame(read_excel("/media/4tb1/Paula/33Br_Nottingham/StatsFusions-Nottingham.xlsx"))
colnames(stats_nottingham)[6] <- "MTT"
stats_nottingham$Fusiones_conf_H[is.na(stats_nottingham$Fusiones_conf_H)] <- 0

boxplots_TFB_MTT(stats = stats_nottingham, group = NA, cohorte = "Nottingham")

fusions_nottingham <- as.data.frame(read_excel("/media/4tb1/Paula/33Br_Nottingham/Todos-FusionReports.xlsx"))
length(which(grepl("kinase", fusions_nottingham$retained_protein_domains))) / nrow(fusions_nottingham) *100
unique(fusions_nottingham$confidence)
length(which(grepl("kinase", fusions_nottingham$retained_protein_domains[fusions_nottingham$confidence == "high"]))) / nrow(fusions_mama61[fusions_mama61$confidence == "high",]) *100

#Breast Normal - SRA ------------------------------------------------------
out <- fusionStats(patients_dir = "/media/4tb2/Daniela/Environ/FusionesActualizado/SRA-NormalBreast/MuestrasSRA",
                   Metadata = NA,
                   group = NA,
                   cohorte = "NormalesBreastSRA")
stats_breast_normal <- read_excel("/media/4tb2/Daniela/Environ/FusionesActualizado/SRA-NormalBreast/MuestrasSRA/StatsFusions_NormalesBreastSRA.xlsx")
boxplots_TFB_MTT(stats = stats_breast_normal, group = NA, cohorte = "NormalesBreastSRA")

fusions_normal <- as.data.frame(read_excel("/media/4tb2/Daniela/Environ/FusionesActualizado/SRA-NormalBreast/MuestrasSRA/Todos-FusionReports_NormalesBreastSRA.xlsx"))
length(which(grepl("kinase", fusions_normal$retained_protein_domains))) / nrow(fusions_normal) *100
unique(fusions_nottingham$confidence)
length(which(grepl("kinase", fusions_nottingham$retained_protein_domains[fusions_nottingham$confidence == "high"]))) / nrow(fusions_mama61[fusions_mama61$confidence == "high",]) *100


# CHINA - SRA ------------------------------------
metadata_china <- Environ_Paula_metadata_global_2025[which(Environ_Paula_metadata_global_2025$categoria_1 == "china"),]
colnames(metadata_china)

  #Esto fue sacado de la documentacion oficial:
meta_china <- read.delim("/media/4tb2/Daniela/Environ/FusionesActualizado/China-SRA/Metadata-China-SraRunTable.txt", sep = ",", stringsAsFactors = FALSE)
colnames(meta_china)[1] <- "ID"
colnames(meta_china)[colnames(meta_china) == "metastasis"] <- "MTT"
meta_china$ID[!meta_china$ID %in% metadata_china$ID]
all(metadata_china$ID %in% meta_china$ID)

metadata_china_completa <- merge(metadata_china[, -which(colnames(metadata_china) == "MTT")], meta_china[, c("ID",  "MTT", "tissue")], by = "ID", all = TRUE)
colnames(metadata_china_completa)
metadata_china_completa$MTT <- ifelse(metadata_china_completa$MTT == "yes", "MET+", "MET-")
out_china <- fusionStats(patients_dir = "/media/4tb2/Daniela/Environ/FusionesActualizado/China-SRA",
                         Metadata = metadata_china_completa,
                         group = "tissue",
                         sobrevida = FALSE)
fusions_china <- out_china[[1]]
stats_china <- out_china[[2]]
genes_china <- out_china[[3]]
TFB_china <- fusions_china[fusions_china$confidence == "high",]

write.xlsx(TFB_china, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_TFB.xlsx")
write.xlsx(metadata_china_completa, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_Metadata.xlsx")
write.xlsx(genes_china, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_Genes.xlsx")
write.xlsx(stats_china, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Breast_China_FusionStats.xlsx")

length(which(grepl("kinase", fusions_china$retained_protein_domains))) / nrow(fusions_china) *100
unique(fusions_china$confidence)
length(which(grepl("kinase", fusions_china$retained_protein_domains[fusions_china$confidence == "high"]))) / nrow(fusions_mama61[fusions_mama61$confidence == "high",]) *100


table(stats_china$MTT)
table(stats_china$Grupo)
library(dplyr)

stats_china$Cohorte <- "China"
stats_breast_normal$MTT <- "-"
stats_breast_normal$Grupo <- "normal breast tissue"

stats_china_normal <- rbind(stats_china[, c("ID", "Fusiones_conf_H", "MTT", "Grupo" ,"Cohorte")],
                            stats_breast_normal[, c("ID", "Fusiones_conf_H", "MTT", "Grupo","Cohorte")])

conteo_mtt_por_grupo <- stats_china %>%
  group_by(Grupo, MTT) %>%
  summarise(n = n(), .groups = "drop")

ggplot(stats_china_normal, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "TFB China + Normales - Breast",
       x = "MTT",
       y = "TFB") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))+
  facet_grid( ~ Grupo)

stats_china %>%
  group_by(Grupo) %>%
  summarise(p_value = wilcox.test(Fusiones_conf_H ~ MTT)$p.value)


# Graficar todas las cohortes juntas: -------------------------------------------

stats_mama61$Cohorte <- "Mama61"
stats_CGD$Cohorte <- "CDGenomics"
stats_tcl$Cohorte <- "TCL"
stats_TCGA$Cohorte <- "TCGA"
stats_nottingham$Cohorte <- "Nottingham"
stats_china$Cohorte <- "China"

colnames(stats_mama61)
stats_conjunto <- rbind(stats_mama61[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_CGD[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_tcl[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_china[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_TCGA[, c("ID", "Fusiones_conf_H", "MTT", "Cohorte")],
                        stats_nottingham[,c("ID", "Fusiones_conf_H", "MTT", "Cohorte")])

any(is.na(stats_conjunto))
write.xlsx(stats_conjunto, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Stats_TodasFusiones.xlsx")


metadata_mama61$Cohorte <- "Mama61"
metadata_CDG$Cohorte <- "CDGenomics"
metadata_tcl$Cohorte <- "TCL"
metadata_tcga$Cohorte <- "TCGA"
metadata_nott$Cohorte <- "Nottingham"
metadata_china$Cohorte <- "China"

metadata_conjunto <- bind_rows(metadata_mama61,
                           metadata_CDG,
                           metadata_tcl,
                           metadata_china,
                           metadata_tcga,
                           metadata_nott)
metadata_conjunto <- metadata_conjunto[, colSums(!is.na(metadata_conjunto)) > 0]
write.xlsx(metadata_conjunto, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/Metadata_TodasCohortes.xlsx")


gen_mama61$Cohorte <- "Mama61"
genes_CGD$Cohorte <- "CDGenomics"
gen
metadata_tcl$Cohorte <- "TCL"
metadata_tcga$Cohorte <- "TCGA"
metadata_nott$Cohorte <- "Nottingham"
metadata_china$Cohorte <- "China"

colnames(stats_conjunto)
TodasFusiones_Metadata_Conjunto <- merge(stats_conjunto[, c("ID", "Fusiones_conf_H" )], metadata_conjunto, by = "ID")
write.xlsx(TodasFusiones_Metadata_Conjunto, file = "/media/4tb2/Daniela/Environ/FusionesActualizado/FusionesyMetadata_TodasCohortes.xlsx")
table(TodasFusiones_Metadata_Conjunto$Cohorte)

ggplot(stats_conjunto, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "Distribución de la Cantidad de Fusiones HIGH por MTT",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal() +
  #scale_y_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 70, by = 10)) +
  facet_wrap(~ Cohorte, ncol = 5)

# Aplicando el test de Wilcoxon por cohorte
stats_conjunto %>%
  group_by(Cohorte) %>%
  summarise(p_value = wilcox.test(Fusiones_conf_H ~ MTT)$p.value)

#Ggplot por subtipo:
ggplot(stats_mama61, aes(x = MTT, y = Fusiones_conf_H, fill = MTT)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  labs(title = "Cantidad de Fusiones HIGH por MTT y subtipo - Mama 61",
       x = "MTT",
       y = "Cantidad de Fusiones") +
  theme_minimal() +
  #scale_y_continuous(limits = c(0, NA)) +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 70, by = 10)) +
  facet_grid( ~ Grupo)  # Facet por Cohorte y Subtipo

# Aplicando el test de Wilcoxon por subtipo
stats_mama61 %>%
  group_by(Grupo) %>%
  summarise(p_value = wilcox.test(Fusiones_conf_H ~ MTT)$p.value)




