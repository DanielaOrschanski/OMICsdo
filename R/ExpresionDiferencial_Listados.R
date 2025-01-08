#' @import ComplexHeatmap
#' @import RColorBrewer
#' @import enrichR
#' @import tidyverse
#' @title Expresion Diferencial Listados
#' @description Generates differential expression analysis.
#' Selects differential genes (depending on p value and lfc) and identify which of them are included in each column in df_listado.
#' Creates heatmaps for:
#' - General DEG
#' - Genes listed in df_listado
#' - Genes DEG included in listes.
#' @param b DGEList
#' @return Control DB updated
#' @examples ControlDB <- AddControl(path_dir ="/home/sam/Patients/123/123.fastq")
#' @export

ExpresionDiferencial_Listados <- function(b, df_listado, Grupos, tejido, lfc = 1, p = 0.05, categoria = "MTT", plotMDS = FALSE, plot_heatmaps = TRUE, porListadoGenes = TRUE) {

  print(table(b$samples[[categoria]]))
  if(Grupos == "") {
    b_G <- b
   #Me quedo con los pacientes que son G3 o G4, elimino los que son 0:
  } else if(any(b$samples[, Grupos] == "0")) {
    b_G <- b[, -which(b$samples[, Grupos] == "0")]
  } else {
    b_G <- b
  }


  print(table(b_G$samples[[categoria]]))

  m <- b_G$samples
  #me quedo con los genes que quiero: ESTO NO SE HACE ACA!!!!!!!!!!!!!!!!!!!!!!!
  #b_G <- b_G[which(b_G$genes$SYMBOL %in% genes_seleccionados),]

  group <- factor(b_G$samples[[categoria]])

  #plot MDS --------------------------------------------------
  if (plotMDS == TRUE) {
    Sample_name <- b_G$samples$ID_ENVIRON
    col.group <- factor(b_G$samples[[categoria]])
    library(RColorBrewer)
    levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
    col.group <- as.character(col.group)

    subtiposL<- b_G$samples$Subtipo
    subtiposL <- ifelse(subtiposL == "Luminal HER2", "HER2", subtiposL)
    plotMDS(cpm(b_G$counts, log=TRUE),
            labels=Sample_name,
            #labels = subtiposL,
            col = col.group,
            dim.plot = c(1,2))

    title(main=sprintf("Sample groups:  %s y %s", unique(b_G$samples[[categoria]])[1], unique(b_G$samples[[categoria]])[2] ))
    legend("topright",
           legend = unique(b_G$samples[[categoria]]),
           fill = unique(col.group),
           pch = 1,
           title = "Categoría")

  }


  #Removing genes lowly expressed genes -------------------------
  Dim1_e <- dim(b_G)
  N <- round(table(b_G$samples[[categoria]])[[1]]/2)
  keep.exprs <- filterByExpr(b_G, group = group, large.n= N)

  #Se queda con los genes con conteos >10 en una cant minima de muestras
  b_G <- b_G[keep.exprs, keep.lib.sizes = FALSE]
  Dim2_e <-dim(b_G)

  # Normalization ---------------------------
  b_G <- calcNormFactors(b_G, method = "TMM")

  # Analisis de expresión diferencial ##############################################

  #1. Design matrix
  design <- model.matrix(~1 + b_G$samples[[categoria]])
  #~1 indica que solo estás incluyendo la intersección

  colnames(design) <- gsub("group","", colnames(design))

  #2. Transform RNA-Seq Data Ready for Linear Modelling
  v <- voom(b_G, design, plot = F)

  #decreasing trend between the means and variances resulting from a combination of technical variation in the sequencing experiment and biological variation amongst the replicate samples from different cell populations.

  #3. Fit linear model for each gene given a series of arrays:
  vfit <- lmFit(v, design)
  #Empirical Bayes Statistics for Differential Expression:
  efit <- treat(vfit, lfc = lfc)

  #Para el informe que quieren:
  #Informe_Mama <- data.frame(b_G$genes)
  Informe_Mama <- as.data.frame(v$E)
  Informe_Mama$p_value <- efit$p.value[,2]
  Informe_Mama$p_adj <- p.adjust(efit$p.value[,2], method="fdr")
  Informe_Mama$stat <- efit$t[,2]
  Informe_Mama$log2FC <- efit$coefficients[,2]
  Informe_Mama$stdev <- efit$stdev.unscaled[,2]
  rownames(Informe_Mama) <- b_G$genes[rownames(Informe_Mama), 2]

  #Para cada columna del resto - ordenadas primero las de grupo A y dsp las del B:

  E_mama <- as.data.frame(v$E)
  E_mama$p_valor <- efit$p.value
  rownames(E_mama) <- b_G$genes[rownames(E_mama), 2]
  #E_mama$ENSEMBL <- b_G$genes[rownames(E_mama), 3]
  #E_mama$SYMBOL <- b_G$genes[rownames(E_mama), 2]

  write.xlsx(E_mama, file = sprintf("/media/4tb1/Daniela/Environ/Solicitud/MatrizExpresion_%s_%s_lfc%s.xlsx", tejido, Grupos, lfc), rowNames = TRUE)
  write.xlsx(Informe_Mama, file = sprintf("/media/4tb1/Daniela/Environ/Informe_%s_%s_lfc%s.xlsx", tejido, Grupos, lfc), rowNames = TRUE)


  # Buscamos un p-valor <0.05 entre los genes segun el group
  ps <- efit$p.value
  message(paste("la cantidad de genes diferenciales son: ", sum(efit$p.value[,2] < p)))
  if(sum(efit$p.value[,2] < p) == 0) {
    return(message("No se encuentran genes diferenciales, probar con un valor lfc más bajo"))
  }


  #####################################################################3
  # HEATMAP CON GENES DIFERENCIALES (independiente de listados) ------
  ################################################################3
  E_mama <- as.data.frame(v$E)
  E_mama$p_valor <- efit$p.value[,2]
  E <- E_mama[efit$p.value[, 2] < p,]
  rownames(E) <- b_G$genes[rownames(E), 2]
  genes_diferenciales <- rownames(E)
  E_gd <- E
  E_gd <- E_gd[order(E_gd$p_valor, decreasing = TRUE),]
  str(E_gd)
  M <- E_gd[, -ncol(E_gd)]
  library(ComplexHeatmap)
  Ms <- t(scale(t(M)))

  #En vez de los nombres de las muestras se pone el nivel del group:
  all(colnames(Ms) == b_G$samples[,1])

  if (tejido == "Mama") {
    b_G$samples$Subtipo <- ifelse(b_G$samples$Subtipo == "HER2", "Luminal HER2", b_G$samples$Subtipo)
    unique(b_G$samples$Subtipo)
  }

  tejido_rep <- rep(tejido, nrow(b_G$samples))
  hm_annot <- ComplexHeatmap::Heatmap(
    Ms,
    column_title = "Expresión solo Genes Diferenciales",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    border =1,
    column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
    row_names_gp = gpar(fontsize = 8),
    show_column_names = TRUE,
    top_annotation = HeatmapAnnotation(
      Pareo = b_G$samples[,Grupos],
      Metastasis = b_G$samples$MTT,
      Subtipo = ifelse(tejido_rep == "Mama", b_G$samples$Subtipo, b_G$samples$Grupo.x)
    )
  )

  if(plot_heatmaps == TRUE) {
    print(hm_annot)
  }


  hm_split <- ComplexHeatmap::Heatmap(
    Ms,
    column_title = "Expresión solo Genes Diferenciales por subtipo",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    border =1,
    column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
    row_names_gp = gpar(fontsize = 8),
    column_split = ifelse(tejido_rep == "Mama", b_G$samples$Subtipo, b_G$samples$Grupo.x),
    show_column_names = TRUE,
    top_annotation = HeatmapAnnotation(
      Pareo =b_G$samples[,Grupos],
      Metastasis = b_G$samples$MTT,
      Subtipo = ifelse(tejido_rep == "Mama", b_G$samples$Subtipo, b_G$samples$Grupo.x)
    )
  )

  if(plot_heatmaps == TRUE) {
    print(hm_split)
  }

  #gd <- as.data.frame(M[genes_diferenciales, ncol(M)])
  #rownames(gd) <- rownames(M[genes_diferenciales,])
  gd <- as.data.frame(E_gd$p_valor)
  gd$Genes <- rownames(E_gd)
  write.xlsx(gd, sprintf("/media/4tb1/Daniela/Environ/Solicitud/GenesDiferenciales_%s_%s_p%s.xlsx", tejido, Grupos, p))


  #####################################################################
  #Enriquecimiento de vias ------------------------------------------------
  #######################################################
  #Lo hago solo con los genes diferenciales totales, independiente de listados

  E <- v$E[efit$p.value[, 2] < p,]
  E <- data.frame(E)
  rownames(E) <- b_G$genes[rownames(E), 2]
  Ms <- t(scale(t(E)))

  # Guardo genes up y genes down ----------
  #Me quedo con las 2 ramas principales y guardo los genes de cada una ----
  dend_H <- hclust(dist(Ms))
  plot(dend_H, cex= 0.6)

  cut_dend <- cutree(dend_H, k=2)
  table(cutree(dend_H, k=2))
  UP <- names(which(cut_dend == 1))
  DOWN <- names(which(cut_dend == 2))

  max_length <- max(length(UP), length(DOWN))

  # Rellena los vectores con NA para igualar sus longitudes
  UP <- c(UP, rep(NA, max_length - length(UP)))
  DOWN <- c(DOWN, rep(NA, max_length - length(DOWN)))

  genes_up_down <- data.frame("Genes UP" = UP, "Genes Down" = DOWN)

  write.xlsx(genes_up_down, sprintf("/media/4tb1/Daniela/Environ/Solicitud/genesUPyDOWN_%s_%s.xlsx", Grupos, categoria))
  #writeLines(names(which(cut_dend == 1)), sprintf("/media/4tb1/Daniela/Environ/Solicitud/genesUP_%s_%s.txt", Grupos, categoria))
  #writeLines(names(which(cut_dend == 2)), sprintf("/media/4tb1/Daniela/Environ/Solicitud/genesDOWN_%s_%s.txt", Grupos, categoria))

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

  plotdown <- plotEnrich(
    down_vias,
    #showTerms = 15,
    #numChar = 40,
    y = "Count",
    orderBy = "P.value",
    xlab = NULL,
    ylab = NULL,
    title = sprintf("Vías metabólicas en - expresados - %s", Grupos)
  )
  print(plotdown)

  plotup <- plotEnrich(
    up_vias,
    #showTerms = 15,
    #numChar = 40,
    y = "Count",
    orderBy = "P.value",
    xlab = NULL,
    ylab = NULL,
    title = sprintf("Vías metabólicas en + expresados - %s", Grupos)
  )

  print(plotup)
  #list_vias[[i]] <- vias
  write.xlsx(vias, sprintf("/media/4tb1/Daniela/Environ/Solicitud/Vias_%s_%s.xlsx", tejido, Grupos), overwrite = TRUE, rowNames= TRUE)



  ######################################################################################################
  #A partir de aca grafico con los genes del listado ###################################################
  ##############################################################################
  plots_listados_split_todos <- list()
  plots_listados_split_selec_dif <- list()
  plots_listados_annot_todos <- list()
  plots_listados_annot_selec_dif <- list()

  if (porListadoGenes == TRUE) {

    InfoListados <- data.frame(matrix(ncol= 6, nrow= 13))

    #  i=1
    #  i=7 #colagenos
    #  i=14 #mergex3

    for ( i in 1:ncol(df_listado)) {
      listado <- colnames(df_listado)[i]
      InfoListados[i, "Listado"] <- listado

      genes_seleccionados <- as.character(na.omit(df_listado[,i][[1]]))
      if( ncol(df_listado) == 1) {
        genes_seleccionados <- df_listado[,1]
      }

      genes_totales <- length(genes_seleccionados)
      InfoListados[i, "Genes_Enlistados"] <- genes_totales

      genes_no_encontrados <- genes_seleccionados[which(!(genes_seleccionados %in% b$genes$SYMBOL))]
      InfoListados[i, "Cant_Genes_NO_Encontrados"] <- ifelse(length(nchar(genes_no_encontrados)) == 0, 0, length(nchar(genes_no_encontrados)))

      genes_filtrados <- genes_seleccionados[which(!(genes_seleccionados %in% b_G$genes$SYMBOL))]
      genes_concatenados <- paste(genes_filtrados, collapse = ", ")
      InfoListados[i, "Genes_Filtrados"] <- genes_concatenados

      #counts_genes_out <- b$counts[which(b$genes$SYMBOL %in% genes_out),]

      #Genes diferenciales
      InfoListados[i, "Cant_Genes_Diferenciales"] <- length(which(genes_seleccionados %in% genes_diferenciales))
      InfoListados[i, "Genes_Diferenciales"] <- paste(genes_seleccionados[which(genes_seleccionados %in% genes_diferenciales)], collapse= ",")

      # Heatmaps con TODOS los genes seleccionados del listado ---------------------
      M <- v$E
      M <- data.frame(M)
      rownames(M) <- b_G$genes[rownames(M), 2]
      M <- M[which(rownames(M) %in% genes_seleccionados),]
      Ms <- t(scale(t(M)))

      all(colnames(Ms) == b_G$samples$ID_ENVIRON)

      if( tejido == "Mama") {
        b_G$samples$Subtipo <- ifelse(b_G$samples$Subtipo == "HER2", "Luminal HER2", b_G$samples$Subtipo)
        unique(b_G$samples$Subtipo)
      }


      hm_annot <- ComplexHeatmap::Heatmap(
        Ms,
        column_title = sprintf("Expresión todos los genes - %s", listado),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        border =1,
        column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
        row_names_gp = gpar(fontsize = 8),
        show_column_names = TRUE,
        top_annotation = HeatmapAnnotation(
          Pareo = b_G$samples[,Grupos],
          Metastasis = b_G$samples$MTT,
          Subtipo = ifelse(tejido_rep == "Mama", b_G$samples$Subtipo, b_G$samples$Grupo.x)
          #df = b_G$samples[,15:ncol(b_G$samples)]
        )
      )
      print(hm_annot)

      hm_split <- ComplexHeatmap::Heatmap(
        Ms,
        column_title = sprintf("Expresión todos los genes por Subtipo - %s", listado),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        border =1,
        column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
        row_names_gp = gpar(fontsize = 8),
        column_split = ifelse(tejido_rep == "Mama", b_G$samples$Subtipo, b_G$samples$Grupo.x),
        show_column_names = TRUE,
        top_annotation = HeatmapAnnotation(
          Pareo = b_G$samples[, Grupos],
          Metastasis = b_G$samples$MTT,
          Subtipo = ifelse(tejido_rep == "Mama", b_G$samples$Subtipo, b_G$samples$Grupo.x)
          #df = b_G$samples[,15:ncol(b_G$samples)]
        )
      )
      print(hm_split)

      plots_listados_annot_todos[[i]] <- hm_annot
      plots_listados_split_todos[[i]] <- hm_split

      #DIFERENCIALES  ---------------------------------------------------------------------
      #M ya tiene solo los seleccionados
      #Veo cuales de los seleccionados estan en los diferenciales
      M <- M[which(rownames(M) %in% genes_diferenciales),]

      library(ComplexHeatmap)
      Ms <- t(scale(t(M)))

      #En vez de los nombres de las muestras se pone el nivel del group:
      all(colnames(Ms) == b_G$samples$ID_ENVIRON)

      if(tejido == "Mama") {
        b_G$samples$Subtipo <- ifelse(b_G$samples$Subtipo == "HER2", "Luminal HER2", b_G$samples$Subtipo)
        unique(b_G$samples$Subtipo)
      }


      hm_annot <- ComplexHeatmap::Heatmap(
        Ms,
        column_title = sprintf("Expresión Genes Diferenciales - %s", listado),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        border =1,
        column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
        row_names_gp = gpar(fontsize = 8),
        show_column_names = TRUE,
        top_annotation = HeatmapAnnotation(
          Pareo = b_G$samples[,Grupos],
          Metastasis = b_G$samples$MTT,
          Subtipo = b_G$samples$Subtipo
        )
      )
      print(hm_annot)

      hm_split <- ComplexHeatmap::Heatmap(
        Ms,
        column_title = sprintf("Expresión Genes Diferenciales por Subtipo - %s", listado),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        border =1,
        column_names_gp = gpar(fontsize = 10),# Ajustar el tamaño de la fuente de los nombres de las columnas
        row_names_gp = gpar(fontsize = 8),
        column_split = b_G$samples$Subtipo,
        show_column_names = TRUE,
        top_annotation = HeatmapAnnotation(
          Pareo =b_G$samples[,Grupos],
          Metastasis = b_G$samples$MTT,
          Subtipo = b_G$samples$Subtipo
        )
      )
      print(hm_split)

      plots_listados_annot_selec_dif[[i]] <- hm_annot
      plots_listados_split_selec_dif[[i]] <- hm_split

      i <- i+1

    }

    InfoListados <- InfoListados[, 7:12]
    write.xlsx(InfoListados, sprintf("/media/4tb1/Daniela/Environ/Solicitud/InfoListados_%s.xlsx", Grupos))

  } else {
    InfoListados <- "No se realizo ningun analisis por listado"
  }

  return(list(gd, vias, E_mama, E_gd, InfoListados, plots_listados_annot_selec_dif, plots_listados_split_selec_dif, plots_listados_annot_todos, plots_listados_split_todos))
}
