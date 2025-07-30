#' @title Run Gene Signatures to MX.
#' @description Evaluates specific genes for proliferation, ca20 and cytolicitc score.
#' @param patients_dir path to the folder that contains one folder for each sample.
#' @param Metadata dataframe that contains information about all the samples. One sample per row.
#' @param group must be a name of one of the columns of the metadata
#' @return Todos_FusionReport and Stats_Fusions: 2 dataframes.
#' @export
#' @import openxlsx
#' @import readxl
#' @import ggpubr
#' @examples fusionStats(patients_dir, Metadata, group = "group")

Correlacion_Prolif <- function(cohorte, TPMt) {

  Prolif_Signature <- read_excel("/media/16TBDisk/Daniela/OMICsdo/data/Prolif_Signature.xlsx")

  if(cohorte == "TrainTPM") {
    mx_train_met_neg_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_train_met_neg_tpm.rds")
    mx_train_met_pos_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_train_met_pos_tpm.rds")
    TPMt <- rbind(mx_train_met_neg_tpm, mx_train_met_pos_tpm)

  } else if (cohorte == "TCLTPM") {
    mx_tcl_met_neg_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_tcl_met_neg_tpm.rds")
    mx_tcl_met_pos_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_tcl_met_pos_tpm.rds")
    TPMt <- rbind(mx_tcl_met_neg_tpm, mx_tcl_met_pos_tpm)

  } else if (cohorte == "CDTPM") {
    mx_cd_met_neg_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_cd_met_neg_tpm.rds")
    mx_cd_met_pos_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_cd_met_pos_tpm.rds")
    TPMt <- rbind(mx_cd_met_neg_tpm, mx_cd_met_pos_tpm)

  } else if (cohorte == "ChinaTPM") {
    TPMt <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_china_todas_TPM_arm.rds")

  } else if (cohorte == "NottTPM") {
    mx_nott_met_neg_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_nott_met_neg_tpm.rds")
    mx_nott_met_pos_tpm <- readRDS("/media/respaldo8t/Daniela/MIXTURE para Environ/mx_val_nott_met_pos_tpm.rds")
    TPMt <- rbind(mx_nott_met_neg_tpm, mx_nott_met_pos_tpm)

  }

  TPMlog <- TPMt+1
  TPMlog <- log(TPMlog)
  prolif <- Prolif_Hoja_1_1_$ZNHIT2
  cols <- match(prolif, colnames(TPMlog))
  cols <- cols[-which(is.na(cols)==T)]
  Prolif_logC <- TPMlog[,cols]

  Prolif_logC <- exp(Prolif_logC)
  Prolif_logC <- as.data.frame(Prolif_logC)

  library(dplyr)
  Prolif_logC <- Prolif_logC %>%
    mutate(Total = rowSums(.))
  N <- nrow(Prolif_logC)
  Prolif_logC$Prolif <- Prolif_logC$Total/N

  #ca20
  t_list <- list(TPMlog)
  rlist <- list()
  ca20<- c("AURKA", "CCNA2", "CCND1", "CCNE2", "CDK1", "CEP63", "CEP152", "E2F1", "E2F2", "LMO4", "MDM2", "MYCN", "NDRG1", "NEK2", "PIN1", "PLK1", "PLK4", "SASS6", "STIL", "TUBG1")

  for (i in 1:1) {
    aux_mx <- t_list[[i]]
    aux_mx <- t(aux_mx)
    aux.expr<- aux_mx-median(aux_mx, na.rm=TRUE)
    aux.expr<- aux.expr/sd(aux.expr, na.rm=TRUE)
    rows <- match(ca20, rownames(aux.expr))
    ca20_aux <- aux.expr[rows,]
    aux.scores<- colSums(ca20_aux, na.rm=TRUE)
    rlist[[i]] <- aux.scores
  }

  ca20_res <- unlist(rlist)

  #cyt score
  cols <- match(c("GZMA", "PRF1"), colnames(TPMt))
  TPM2 <- TPMt[,cols]
  cyt_res <- exp((log(TPM2[,1] + 1)+log(TPM2[,2]+1))/2)

  all(names(cyt_res) == names(ca20_res))
  all(names(cyt_res) == rownames(Prolif_logC))

  Resultados <- cbind("ID" = names(cyt_res), cyt_res, ca20_res, Prolif_logC$Prolif)
  colnames(Resultados)[4] <- "prolif"
  colnames(metadata_completa)
  Resultados <- merge(Resultados, metadata_completa[, c("ID", "Fusiones_conf_H", "MTT")], by ="ID")
  colnames(Resultados)[5] <- "Fidx"

  #Graficos:
  library(tidyr)
  library(ggplot2)

  Resultados_long <- pivot_longer(Resultados, cols = 2:4,
                                  names_to = "Variable", values_to = "Valor")
  Resultados_long$Valor <- as.numeric(as.character(Resultados_long$Valor))

  # Calcular correlaciones
  correlaciones <- Resultados_long %>%
    group_by(MTT, Variable) %>%
    summarise(
      correlacion = cor(Fidx, Valor, use = "complete.obs"),
      .groups = "drop"
    ) %>%
    mutate(correlacion = round(correlacion, 3),
           label = paste0(MTT, ": r = ", correlacion))

  rangos <- Resultados_long %>%
    group_by(Variable) %>%
    summarise(
      x_max = mean(Valor, na.rm = TRUE),
      y_max = max(Fidx, na.rm = TRUE),
      y_min = min(Fidx, na.rm = TRUE),
      .groups = "drop"
    )

  # Unimos rangos con correlaciones para definir posición
  cor_pos <- correlaciones %>%
    left_join(rangos, by = "Variable") %>%
    mutate(x_pos = x_max ,
           # Para separar los dos grupos en y, usamos y_max para MET+ y y_min para MET-
           y_pos = ifelse(MTT == "MET+", y_max+5, y_max+5 + 0.15 * (y_max - y_min))
    )

  g <- ggplot(Resultados_long, aes(x = Valor, y = Fidx, color = MTT)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(group = MTT, color = MTT)) +
    geom_text(data = cor_pos, aes(x = x_pos, y = y_pos, label = label, color = MTT),
              inherit.aes = FALSE, hjust = 0, size = 4, fontface = "bold") +
    facet_wrap(~ Variable, scales = "free_x", nrow = 1) +
    scale_y_continuous(labels = NULL)+
    theme_minimal() +
    labs(x = "Valor de Variable", y = "Fidx", color = "Grupo MTT",
         title = sprintf("Fidx por MTT - %s", cohorte)) +
    theme(strip.text = element_text(size = 12, face = "bold"),
          panel.clip = "on")

  print(g)

}
