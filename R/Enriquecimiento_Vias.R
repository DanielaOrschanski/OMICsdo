#' @title Enriquecimiento vias
#' @description It gein MitoR DataBase.
#' @param b_G Either the BAMfastq files
#' @param categoria Filters applied to the anded by GATK
#' @param Grupos visuf
#' @param p pvalue
#' @param lfc log fold change
#' @param genes_seleccionados list of kfernkdf name of the list
#' @param listado ofsijkkds
#' @return vias egvfnjrnfk
#' @export
#' @examples genes_seleccionados <- Solicitud_2Br_DEG_en_G3vsG4$mergex3
#' vias <- Enriquecimiento_Vias(b_G = b_G, Grupos = "G18G19", genes_seleccionados = genes_seleccionados, listado = "mergex3")

Enriquecimiento_Vias <- function(b_G, categoria = "MTT", Grupos, p = 0.05, lfc = 1, genes_seleccionados = NULL, listado = "") {

  #1. Design matrix
  design <- model.matrix(~1 + b_G$samples[[categoria]])
  #~1 indica que solo estás incluyendo la intersección

  colnames(design) <- gsub("group","", colnames(design))

  #2. Transform RNA-Seq Data Ready for Linear Modelling
  v <- voom(b_G, design, plot = F)

  #decreasing trend between the means and variances resulting from a combination of technical variation in the sequencing experiment and biological variation amongst the replicate samples from different cell populations.

  #3. Fit linear model for each gene given a series of arrays:
  vfit <- lmFit(v, design)
  efit <- treat(vfit, lfc = lfc)

  E_mama <- as.data.frame(v$E)
  E_mama$p_valor <- efit$p.value
  rownames(E_mama) <- b_G$genes[rownames(E_mama), 2]

  ps <- efit$p.value
  message(paste("la cantidad de genes diferenciales son: ", sum(efit$p.value[,2] < p)))

  E <- v$E[efit$p.value[, 2] < p,]
  E <- data.frame(E)
  rownames(E) <- b_G$genes[rownames(E), 2]

  if(!(is.null(genes_seleccionados))) {
    E <- E[which(rownames(E) %in% genes_seleccionados),]
  }

  Ms <- t(scale(t(E)))

  # Guardo genes up y genes down ----------
  #Me quedo con las 2 ramas principales y guardo los genes de cada una ----
  dend_H <- hclust(dist(Ms))
  plot(dend_H, cex= 0.6)

  cut_dend <- cutree(dend_H, k=2)
  table(cutree(dend_H, k=2))
  UP <- names(which(cut_dend == 1))
  DOWN <- names(which(cut_dend == 2))

  #writeLines(names(which(cut_dend == 1)), sprintf("/media/4tb1/Daniela/Environ/Solicitud/genesUP_%s_%s_%s.txt", Grupos, categoria, listado))
  #writeLines(names(which(cut_dend == 2)), sprintf("/media/4tb1/Daniela/Environ/Solicitud/genesDOWN_%s_%s_%s.txt", Grupos, categoria, listado))

  #genes_UP[[i]] <- UP
  #genes_DOWN[[i]] <- DOWN

  #Vinculación de genes con vías metabolicas ------------

  library(enrichR)
  #install.packages("enrichR")
  dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human")

  library(tidyverse)

  enriched_up <- enrichr(UP, dbs)
  enriched_down <- enrichr(DOWN, dbs)

  #Genero un solo dataframe con toda la info de todas las db ------
  down_vias <- do.call("rbind", enriched_down)
  down_vias <- down_vias[which(down_vias$P.value < 0.05),]
  down_vias <- down_vias[order(down_vias$P.value),]

  up_vias <- do.call("rbind", enriched_up)
  up_vias <- up_vias[which(up_vias$P.value < 0.05),]
  up_vias <- up_vias[order(up_vias$P.value),]

  up_vias$type <- "Up_regulated"
  down_vias$type <- "Down_regulated"
  vias <- rbind(down_vias, up_vias)
  vias$cohort <- "Breast"
  vias$Pareo <- Grupos

  plotEnrich(
    down_vias,
    #showTerms = 15,
    #numChar = 40,
    y = "Count",
    orderBy = "P.value",
    xlab = NULL,
    ylab = NULL,
    title = "Vías metabólicas en - expresados"
  )

  plotEnrich(
    up_vias,
    #showTerms = 15,
    #numChar = 40,
    y = "Count",
    orderBy = "P.value",
    xlab = NULL,
    ylab = NULL,
    title = "Vías metabólicas en + expresados"
  )

  #list_vias[[i]] <- vias
  write.xlsx(vias, sprintf("/media/4tb1/Daniela/Environ/Solicitud/Vias_%s_%s_%s.xlsx", Grupos, listado, categoria), overwrite = TRUE, rowNames= TRUE)
  return(vias)

}
