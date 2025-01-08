
cox_porsubtipo <- function(mx, Grupos, listado = "", porSubtipo = TRUE) {

  if(porSubtipo == TRUE) {
    subtipos <- unique(mx$subtipo)
  } else {
    subtipos <- unique(mx$subtipo)[1]
  }


  for (subtipo in subtipos) {
    #subtipo <- subtipos[2]
    #subtipo <- subtipos
    print(subtipo)

    #Para que se haga tanto por subtipo como independiente de subtipo
    if(porSubtipo == TRUE) {
      message("El analsis se va a realizar POR subtipos")
      mx_sub <- mx[which(mx$subtipo == subtipo),]
    } else {
      message("El analisis se va a realizar INDEPENDIENTE de subtipos")
      mx_sub <- mx
      subtipo <- ""
    }


    #Bootstrap a 0 -------------------------
    test_0_bootstrap <- bootstrap_0_lasso(mx_sub)
    bootstrap_0_stats <- test_0_bootstrap[[1]]
    bootstrap_0_models <- test_0_bootstrap[[2]]

    bootstrap_0_genes <- test_0_bootstrap[[3]]
    bootstrap_0_genes[bootstrap_0_genes != 0] <- 1
    bootstrap_0_genes <- cbind(bootstrap_0_genes, BTT0 = rowSums(bootstrap_0_genes))
    summary(colSums(bootstrap_0_genes[,-101]))

    bootstrap_0_genes <- bootstrap_0_genes[order(bootstrap_0_genes[, 101], decreasing = T),]

    #Bootstrap median -------------------------
    test_bootstrap <- bootstrap_median_lasso(mx_sub)
    bootstrap_stats <- test_bootstrap[[1]]
    bootstrap_models <- test_bootstrap[[2]]

    bootstrap_genes <- test_bootstrap[[3]]
    bootstrap_genes[bootstrap_genes != 0] <- 1
    bootstrap_genes <- cbind(bootstrap_genes, BTTM = rowSums(bootstrap_genes))
    summary(colSums(bootstrap_genes[,-101]))

    bootstrap_genes <- bootstrap_genes[order(bootstrap_genes[, 101], decreasing = T),]

    #Custom0 ------------------------------
    #if(custom == TRUE) {
    test_iter_0 <- iter_lasso_0(mx_sub)

    test_iter_0_stats <- test_iter_0[[1]]
    test_iter_0_models <- test_iter_0[[2]]

    glmnet_0_genes <- test_iter_0[[3]]
    glmnet_0_genes[glmnet_0_genes != 0] <- 1
    glmnet_0_genes <- cbind(glmnet_0_genes, custom0 = rowSums(glmnet_0_genes))
    summary(colSums(glmnet_0_genes[,-101]))
    glmnet_0_genes <- glmnet_0_genes[order(glmnet_0_genes[,101], decreasing=T),]


    #Custom median ------------------------------
    test_iter <- iter_lasso_median(mx_sub)

    test_iter_stats <- test_iter[[1]]
    test_iter_models <- test_iter[[2]]
    glmnet_genes <- test_iter[[3]]
    glmnet_genes[glmnet_genes != 0] <- 1
    glmnet_genes <- cbind(glmnet_genes, custom0 = rowSums(glmnet_genes))
    summary(colSums(glmnet_genes[,-101]))

    glmnet_genes <- glmnet_genes[order(glmnet_genes[,101], decreasing=T),]

    #}

    #LOO 0 -------------------------------------

    test_LOO0 <- LOO_0_lasso(mx_sub)

    LOO0_stats <- test_LOO0[[1]]
    LOO0_models <- test_LOO0[[2]]

    LOO0_genes <- test_LOO0[[3]]
    LOO0_genes[LOO0_genes != 0] <- 1
    LOO0_genes <- cbind(LOO0_genes, LOO0 = rowSums(LOO0_genes))
    summary(colSums(LOO0_genes[,-101]))

    LOO0_genes <- LOO0_genes[order(LOO0_genes[,101], decreasing=T),]

    #LOO median -------------------------------------

    test_LOOM <- LOO_median_lasso(mx_sub)

    LOOM_stats <- test_LOOM[[1]]
    LOOM_models <- test_LOOM[[2]]

    LOOM_genes <- test_LOOM[[3]]
    LOOM_genes[LOOM_genes != 0] <- 1
    LOOM_genes <- cbind(LOOM_genes, LOOM = rowSums(LOOM_genes))
    summary(colSums(LOOM_genes[,-101]))

    LOOM_genes <- LOOM_genes[order(LOOM_genes[,101], decreasing=T),]



    ###################################################3
    LOOM <- LOOM_genes[order(rownames(LOOM_genes)), ]
    LOO0 <- LOO0_genes[order(rownames(LOO0_genes)), ]

    custom0 <- glmnet_0_genes[order(rownames(glmnet_0_genes)), ]
    customM <- glmnet_genes[order(rownames(glmnet_genes)), ]

    BOOT0 <- bootstrap_0_genes[order(rownames(bootstrap_0_genes)), ]
    BOOTM <- bootstrap_genes[order(rownames(bootstrap_genes)), ]


    # Extrae la columna 101 de cada dataframe
    LOOM <- LOOM[, 101, drop = FALSE]
    LOO0 <- LOO0[, 101, drop = FALSE]

    custom0 <- custom0[, 101, drop = FALSE]
    customM <- customM[, 101, drop = FALSE]


    BOOT0 <- BOOT0[, 101, drop = FALSE]
    BOOTM <- BOOTM[, 101, drop = FALSE]

    # Combina las columnas extraídas en un nuevo dataframe
    combined_df <- cbind(Gen = rownames(BOOT0), BOOT0, BOOTM, LOO0, LOOM, custom0, customM = customM)

    combined_df <- as.data.frame(combined_df)
    combined_df[,-1] <- lapply(combined_df[,-1], as.numeric)
    colnames(combined_df)[6:7] <- c("custom0", "customM")
    combined_df$Promedio <- round(rowMeans(combined_df[, -1], na.rm = TRUE), 2)

    library(openxlsx)
    write.xlsx(combined_df, sprintf("~/Daniela/Environ/CoxLasso_GeneFreq-%s-%s-%s.xlsx", Grupos, subtipo, listado))
  }
}

###################################################################
# VENN COMPARANDO LOS TOP GENES ELEGIDOS EN CADA CASO #############
###################################################################

#cox_G18G19 <- list(CoxLasso_GeneFreq_G18G19_GenesDiferenciales, CoxLasso_GeneFreq_G18G19_HER2_,
#                   CoxLasso_GeneFreq_G18G19_Luminal_A_, CoxLasso_GeneFreq_G18G19_Luminal_B_,
#                   CoxLasso_GeneFreq_G18G19_Triple_Negativo_)
#out2 <- findIntersect(cox = cox_G18G19, umbral = 1, local_path = "~/Daniela/Environ/G18G19-dif", tecnique = "Promedio", Grupos = "G18G19", listado = "")

findIntersect <- function(cox, local_path, umbral, tecnique = "", Grupos = "", listado = "") {

  filtered_genes <- lapply(cox, function(df) {
    df_filtered  <- df[df[[tecnique]] >= umbral, "Gen"]
    df_filtered <- as.character(unlist(df_filtered))
  })
  names(filtered_genes) <- c("G18G19", "HER2", "Luminal A", "Luminal B", "Triple Negativo")
  #filtered_genes <- filtered_genes[sapply(filtered_genes, length) > 0]

  library(VennDiagram)

  colors <- c("red", "blue", "green", "purple", "orange")
  png(filename = sprintf("%s/venn_%s.png", local_path, tecnique), width = 800, height = 600, res = 100) # Ajustar tamaño y resolución según sea necesario

  venn.plot <- venn.diagram(
    x = filtered_genes,
    category.names = names(filtered_genes),
    filename = NULL,
    output = TRUE,
    col = "transparent",                # Sin bordes
    fill = colors[1:length(filtered_genes)],  # Colores de relleno
    alpha = 0.50,                       # Transparencia del relleno
    cex = 1,                          # Tamaño de la fuente de los números
    cat.cex = 1,
    cat.pos = 0,
    cat.dist = 0.05,
    margin = 0,
    cat.col = colors[1:length(filtered_genes)] # Colores de las etiquetas
  )

  grid.draw(venn.plot)
  dev.off()

  library(grid)
  library(ggplot2)


  # Extraer los nombres de las intersecciones -----------------------------------------
  intersection_list <- list()
  # Calcular intersecciones entre pares
  for (i in 1:(length(filtered_genes) - 1)) {
    for (j in (i + 1):length(filtered_genes)) {
      name <- paste(names(filtered_genes)[i], "&", names(filtered_genes)[j])
      intersection_list[[name]] <- intersect(filtered_genes[[i]], filtered_genes[[j]])
    }
  }

  # Calcular intersecciones entre tríos
  for (i in 1:(length(filtered_genes) - 2)) {
    for (j in (i + 1):(length(filtered_genes) - 1)) {
      for (k in (j + 1):length(filtered_genes)) {
        name <- paste(names(filtered_genes)[i], "&", names(filtered_genes)[j], "&", names(filtered_genes)[k])
        intersection_list[[name]] <- Reduce(intersect, list(filtered_genes[[i]], filtered_genes[[j]], filtered_genes[[k]]))
      }
    }
  }

  # Calcular la intersección de todos los conjuntos
  all_intersections_name <- paste(names(filtered_genes), collapse = " & ")
  intersection_list[[all_intersections_name]] <- Reduce(intersect, filtered_genes)

  # Convertir la lista de intersecciones en un data.frame
  intersection_df <- data.frame(
    Intersection = names(intersection_list),
    Genes = sapply(intersection_list, function(genes) {
      paste(genes, collapse = ", ")  # Concatenar genes en una sola cadena separada por comas
    }),
    stringsAsFactors = FALSE
  )

  library(openxlsx)
  write.xlsx(intersection_df, file = sprintf("%s/Interseccion-Genes_%s_%s_por%s%s.xlsx", local_path, Grupos, listado, tecnique, umbral))
  return(intersection_df)
}

#################################################################
# FUNCIONES DE COX Y LASSO - DISTINTAS TECNICAS DE SELECCION DE MUESTRAS
################################################################################

bootstrap_0_lasso <- function(mx) {
  #stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  #colnames(stat_iter) <- c("specificity", "sensitivity", "accuracy", "ppv", "npv", "tn", "tp", "fn", "fp","D","F1")

  #yo:
  stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  colnames(stat_iter) <- c( "Sensitivity", "Specificity",  "ppv", "npv", "Precision", "Recall", "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")

  list_model <- list()
  list_genes <- list()

  data <- mx
  expresion <- as.matrix(data[, -c(1:3)])
  metadata <- as.matrix(data[,c(1:3)])

  tiempo <- data$tiempo
  evento <- data$evento

  pos_pat <- rownames(data[which(data$evento == 1),])
  neg_pat <- rownames(data[which(data$evento == 0),])

  for(x in 1:100) {
    #Saca con reemplazo la misma cantidad de pacientes pos que hay:
    id_pos <- sample(pos_pat, replace = TRUE)
    id_neg <- sample(neg_pat, replace = TRUE)
    #train <- as.numeric(c(id_pos, id_neg))
    train <- c(id_pos, id_neg)

    expresion_train <- expresion[train,]
    #tiempo_train <- data$tiempo[train]
    tiempo_train <- data[train, "tiempo"]
    #evento_train <- data$evento[train]
    evento_train <- data[train, "evento"]

    col_test <- unique(train)

    expresion <- as.data.frame(expresion)
    expresion_test <- expresion[-which(rownames(expresion) %in% col_test),]

    length(col_test) + nrow(expresion_test)

    #install.packages("verification")

    surv_obj <- Surv(tiempo_train, evento_train)
    cox_lasso <- cv.glmnet(as.matrix(expresion_train), surv_obj, family = "cox", alpha = 1)

    #Hacer esto pero apra un clasificador

    best_lam <- cox_lasso$lambda.min
    lasso_best <- glmnet(expresion_train, surv_obj, family = "cox", alpha = 1, lambda = best_lam)
    #hacer lasso_best pero con lambda infinito

    model_matrix <- as.matrix(lasso_best[["beta"]])

    #riesgo_lineal <- predict(lasso_best, newx = expresion_test, type = "link")
    riesgo_lineal <- predict(lasso_best, newx = as.matrix(expresion_test), type = "link")

    grupo_riesgo <- ifelse(riesgo_lineal > 0, "High Risk", "Low Risk")

    #agregado por mi:
    #hay que hacer algo si es una sola:
    if(nrow(grupo_riesgo) == 1) {
      meta <- metadata[-which(rownames(metadata) %in% train), c("subtipo","tiempo", "evento")]
      id <- rownames(metadata)[!(rownames(metadata) %in% train)]
      g <- as.data.frame(t(meta))
      rownames(g) <- id
      g$s0
      grupo_riesgo <- cbind(g, grupo_riesgo)
    } else {
      grupo_riesgo <- cbind(metadata[-which(rownames(metadata) %in% train), c("subtipo","tiempo", "evento")], grupo_riesgo)
    }

    grupo_riesgo <- as.data.frame(grupo_riesgo)
    grupo_riesgo$predicho <- ifelse(grupo_riesgo$s0 == "High Risk", 1, 0)

    #Lo que hace guada:
    #tab <- table(grupo_riesgo$evento, grupo_riesgo$s0)
    #stats_tab <- table.stats(tab)
    #stat_iter[x,] <- stats_tab

    #lo que cambie yo:
    library(caret)
    print(x)
    grupo_riesgo$predicho <- factor(grupo_riesgo$predicho, levels = c("0", "1"))
    grupo_riesgo$evento <- factor(grupo_riesgo$evento, levels = c("0", "1"))

    cf0 <- confusionMatrix(as.factor(grupo_riesgo$predicho), as.factor(grupo_riesgo$evento), mode = "everything", positive = "1")
    stat_iter[x,] <- cf0$byClass

    list_model[[x]] <- grupo_riesgo
    list_genes[[x]] <- model_matrix
  }

  list_genes <- do.call(cbind, list_genes)
  colnames(list_genes) <- 1:100

  return(list(stat_iter, list_model, list_genes))
}


bootstrap_median_lasso <- function(mx) {
  #stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  #colnames(stat_iter) <- c("specificity", "sensitivity", "accuracy", "ppv", "npv", "tn", "tp", "fn", "fp","D","F1")

  #yo:
  stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  colnames(stat_iter) <- c( "Sensitivity", "Specificity",  "ppv", "npv", "Precision", "Recall", "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")

  list_model <- list()
  list_genes <- list()

  data <- mx

  expresion <- as.matrix(data[, -c(1:3)])
  metadata <- as.matrix(data[,c(1:3)])

  tiempo <- data$tiempo
  evento <- data$evento

  pos_pat <- rownames(data[which(data$evento == 1),])
  neg_pat <- rownames(data[which(data$evento == 0),])

  for(x in 1:100) {
    #Saca con reemplazo la misma cantidad de pacientes pos que hay:
    id_pos <- sample(pos_pat, replace = TRUE)
    id_neg <- sample(neg_pat, replace = TRUE)
    #train <- as.numeric(c(id_pos, id_neg))
    train <- c(id_pos, id_neg)

    expresion_train <- expresion[train,]
    #tiempo_train <- data$tiempo[train]
    tiempo_train <- data[train, "tiempo"]
    #evento_train <- data$evento[train]
    evento_train <- data[train, "evento"]

    col_test <- unique(train)

    expresion <- as.data.frame(expresion)
    expresion_test <- expresion[-which(rownames(expresion) %in% col_test),]

    length(col_test) + nrow(expresion_test)

    #install.packages("verification")
    library(glmnet)
    library(survival)
    library(verification)

    surv_obj <- Surv(tiempo_train, evento_train)
    cox_lasso <- cv.glmnet(as.matrix(expresion_train), surv_obj, family = "cox", alpha = 1)
    best_lam <- cox_lasso$lambda.min
    lasso_best <- glmnet(expresion_train, surv_obj, family = "cox", alpha = 1, lambda = best_lam)

    model_matrix <- as.matrix(lasso_best[["beta"]])
    riesgo_train <- predict(lasso_best, newx = as.matrix(expresion_train), type = "link")
    umbral <- median(riesgo_train)

    #riesgo_lineal <- predict(lasso_best, newx = expresion_test, type = "link")
    riesgo_lineal <- predict(lasso_best, newx = as.matrix(expresion_test), type = "link")
    grupo_riesgo <- ifelse(riesgo_lineal > umbral, "High Risk", "Low Risk")

    #agregado por mi:
    #hay que hacer algo si es una sola:
    if(nrow(grupo_riesgo) == 1) {
      meta <- metadata[-which(rownames(metadata) %in% train), c("subtipo","tiempo", "evento")]
      id <- rownames(metadata)[!(rownames(metadata) %in% train)]
      g <- as.data.frame(t(meta))
      rownames(g) <- id
      g$s0
      grupo_riesgo <- cbind(g, grupo_riesgo)
    } else {
      grupo_riesgo <- cbind(metadata[-which(rownames(metadata) %in% train), c("subtipo","tiempo", "evento")], grupo_riesgo)
    }

    grupo_riesgo <- as.data.frame(grupo_riesgo)
    grupo_riesgo$predicho <- ifelse(grupo_riesgo$s0 == "High Risk", 1, 0)

    #Lo que hace guada:
    #tab <- table(grupo_riesgo$evento, grupo_riesgo$s0)
    #stats_tab <- table.stats(tab)
    #stat_iter[x,] <- stats_tab

    #lo que cambie yo:
    library(caret)
    grupo_riesgo$predicho <- factor(grupo_riesgo$predicho, levels = c("0", "1"))
    grupo_riesgo$evento <- factor(grupo_riesgo$evento, levels = c("0", "1"))

    cf0 <- confusionMatrix(as.factor(grupo_riesgo$predicho), as.factor(grupo_riesgo$evento), mode = "everything", positive = "1")
    stat_iter[x,] <- cf0$byClass

    list_model[[x]] <- grupo_riesgo
    list_genes[[x]] <- model_matrix
  }

  list_genes <- do.call(cbind, list_genes)
  colnames(list_genes) <- 1:100

  return(list(stat_iter, list_model, list_genes))
}


#####################################################################
iter_lasso_median <- function(mx) {
  #stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  #colnames(stat_iter) <- c("specificity", "sensitivity", "accuracy", "ppv", "npv", "tn", "tp", "fn", "fp","D","F1")

  stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  colnames(stat_iter) <- c( "Sensitivity", "Specificity",  "ppv", "npv", "Precision", "Recall", "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")

  list_model <- list()
  list_genes <- list()

  data <- mx
  expresion <- as.matrix(data[, -c(1:3)])
  metadata <- as.matrix(data[,c(1:3)])

  tiempo <- data$tiempo
  evento <- data$evento

  for(x in 1:100) {
    #train = sample(1:nrow(expresion), nrow(expresion)/3*2)
    train = sample(1:nrow(expresion), nrow(expresion)*0.80)
    expresion_train <- expresion[train,]
    tiempo_train <- data$tiempo[train]
    evento_train <- data$evento[train]

    expresion_test <- expresion[-train,]

    surv_obj <- Surv(tiempo_train, evento_train)
    cox_lasso <- cv.glmnet(expresion_train, surv_obj, family = "cox", alpha = 1)
    best_lam <- cox_lasso$lambda.min
    lasso_best <- glmnet(expresion_train, surv_obj, family = "cox", alpha = 1, lambda = best_lam)

    model_matrix <- as.matrix(lasso_best[["beta"]])
    riesgo_train <- predict(lasso_best, newx = as.matrix(expresion_train), type = "link")
    umbral <- median(riesgo_train)

    riesgo_lineal <- predict(lasso_best, newx = expresion_test, type = "link")
    grupo_riesgo <- ifelse(riesgo_lineal > umbral, "High Risk", "Low Risk")
    #grupo_riesgo <- cbind(METADATA_Actualizada[-train, c("ID.ENVIRON", "tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- cbind(metadata[-train, c("subtipo","tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- as.data.frame(grupo_riesgo)

    #grupo_riesgo <- cbind(METADATA_Actualizada[-which(METADATA_Actualizada$ID.ENVIRON %in% train), c("ID.ENVIRON", "tiempo", "evento")], grupo_riesgo)
    grupo_riesgo$predicho <- ifelse(grupo_riesgo$s0 == "High Risk", 1, 0)


    #tab <- table(grupo_riesgo$evento, grupo_riesgo$s0)
    #stats_tab <- table.stats(tab)
    #stat_iter[x,] <- stats_tab

    grupo_riesgo$predicho <- factor(grupo_riesgo$predicho, levels = c("0", "1"))
    grupo_riesgo$evento <- factor(grupo_riesgo$evento, levels = c("0", "1"))

    cf0 <- confusionMatrix(as.factor(grupo_riesgo$predicho), as.factor(grupo_riesgo$evento), mode = "everything", positive = "1")
    stat_iter[x,] <- cf0$byClass

    list_model[[x]] <- grupo_riesgo
    list_genes[[x]] <- model_matrix
  }

  list_genes <- do.call(cbind, list_genes)
  colnames(list_genes) <- 1:100

  return( list(stat_iter, list_model, list_genes))
}

########################################################################
iter_lasso_0 <- function(mx) {
  #stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  #colnames(stat_iter) <- c("specificity", "sensitivity", "accuracy", "ppv", "npv", "tn", "tp", "fn", "fp","D","F1")

  stat_iter <- data.frame(matrix(ncol = 11, nrow = 100))
  colnames(stat_iter) <- c( "Sensitivity", "Specificity",  "ppv", "npv", "Precision", "Recall", "F1", "Prevalence", "Detection_Rate", "Detection_Prevalence", "Balanced_Accuracy")

  list_model <- list()
  list_genes <- list()

  data <- mx
  expresion <- as.matrix(data[, -c(1:3)])
  metadata <- as.matrix(data[,c(1:3)])

  tiempo <- data$tiempo
  evento <- data$evento

  for(x in 1:100) {
    #train = sample(1:nrow(expresion), nrow(expresion)/3*2)
    train = sample(1:nrow(expresion), nrow(expresion)*0.80)
    expresion_train <- expresion[train,]
    tiempo_train <- data$tiempo[train]
    evento_train <- data$evento[train]

    expresion_test <- expresion[-train,]

    surv_obj <- Surv(tiempo_train, evento_train)
    cox_lasso <- cv.glmnet(expresion_train, surv_obj, family = "cox", alpha = 1)
    best_lam <- cox_lasso$lambda.min
    lasso_best <- suppressWarnings(glmnet(expresion_train, surv_obj, family = "cox", alpha = 1, lambda = best_lam))

    model_matrix <- as.matrix(lasso_best[["beta"]])

    riesgo_lineal <- predict(lasso_best, newx = expresion_test, type = "link")
    grupo_riesgo <- ifelse(riesgo_lineal > 0, "High Risk", "Low Risk")
    #grupo_riesgo <- cbind(METADATA_Actualizada[-train, c("ID.ENVIRON", "tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- cbind(metadata[- train, c("subtipo","tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- as.data.frame(grupo_riesgo)

    #grupo_riesgo <- cbind(METADATA_Actualizada[-which(METADATA_Actualizada$ID.ENVIRON %in% train), c("ID.ENVIRON", "tiempo", "evento")], grupo_riesgo)
    grupo_riesgo$predicho <- ifelse(grupo_riesgo$s0 == "High Risk", 1, 0)


    #tab <- table(grupo_riesgo$evento, grupo_riesgo$s0)
    #stats_tab <- table.stats(tab)
    #stat_iter[x,] <- stats_tab
    grupo_riesgo$predicho <- factor(grupo_riesgo$predicho, levels = c("0", "1"))
    grupo_riesgo$evento <- factor(grupo_riesgo$evento, levels = c("0", "1"))

    cf0 <- confusionMatrix(as.factor(grupo_riesgo$predicho), as.factor(grupo_riesgo$evento), mode = "everything", positive = "1")
    stat_iter[x,] <- cf0$byClass

    list_model[[x]] <- grupo_riesgo
    list_genes[[x]] <- model_matrix
  }

  list_genes <- do.call(cbind, list_genes)
  colnames(list_genes) <- 1:100

  return( list(stat_iter, list_model, list_genes))
}
###############################################################
LOO_0_lasso <- function(mx) {
  stat_iter <- data.frame(matrix(ncol = 4, nrow = 100))
  colnames(stat_iter) <- c("tn", "tp", "fn", "fp")

  list_model <- list()
  list_genes <- list()

  data <- mx
  expresion <- as.matrix(data[, -c(1:3)])
  metadata <- as.matrix(data[,c(1:3)])

  tiempo <- data$tiempo
  evento <- data$evento

  pos_pat <- rownames(data[which(data$evento == 1),])
  neg_pat <- rownames(data[which(data$evento == 0),])

  for(x in 1:100) {
    id_pos <- sample(pos_pat, 1)
    id_neg <- sample(neg_pat, 1)
    #test <- as.numeric(c(id_pos, id_neg))
    test <- c(id_pos, id_neg)

    expresion <- as.data.frame(expresion)
    expresion_train <- expresion[-which(rownames(expresion) %in% test),]
    tiempo_train <- data[-which(rownames(expresion) %in% test), "tiempo"]
    evento_train <- data[-which(rownames(expresion) %in% test), "evento"]

    #tiempo_train <- data$tiempo[-test]
    #evento_train <- data$evento[-test]

    expresion_test <- expresion[test,]

    surv_obj <- Surv(tiempo_train, evento_train)
    cox_lasso <- cv.glmnet(as.matrix(expresion_train), surv_obj, family = "cox", alpha = 1)
    best_lam <- cox_lasso$lambda.min
    lasso_best <- suppressWarnings(glmnet(expresion_train, surv_obj, family = "cox", alpha = 1, lambda = best_lam))

    model_matrix <- as.matrix(lasso_best[["beta"]])

    riesgo_lineal <- predict(lasso_best, newx = as.matrix(expresion_test), type = "link")
    grupo_riesgo <- ifelse(riesgo_lineal > 0, "High Risk", "Low Risk")
    #grupo_riesgo <- cbind(metadata[test,],grupo_riesgo)
    #grupo_riesgo <- cbind(METADATA_Actualizada[which(METADATA_Actualizada$ID.ENVIRON %in% test), c("ID.ENVIRON", "tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- cbind(metadata[which(rownames(metadata) %in% test), c("subtipo","tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- as.data.frame(grupo_riesgo)

    grupo_riesgo$predicho <- ifelse(grupo_riesgo$s0 == "High Risk", 1, 0)


    tab <- table(grupo_riesgo$evento, grupo_riesgo$s0)
    stat_iter$tp[x] <- tab[2]
    stat_iter$tn[x] <- tab[3]
    stat_iter$fn[x] <- tab[4]
    stat_iter$fp[x] <- tab[1]
    list_model[[x]] <- grupo_riesgo
    list_genes[[x]] <- model_matrix
  }

  stat_iter[is.na(stat_iter)] <- 0
  list_genes <- do.call(cbind, list_genes)
  colnames(list_genes) <- 1:100
  result <- list(stat_iter, list_model, list_genes)
  return(result)
}

#############################################################3
LOO_median_lasso <- function(mx) {
  stat_iter <- data.frame(matrix(ncol = 4, nrow = 100))
  colnames(stat_iter) <- c("tn", "tp", "fn", "fp")

  list_model <- list()
  list_genes <- list()

  data <- mx
  expresion <- as.matrix(data[, -c(1:3)])
  metadata <- as.matrix(data[,c(1:3)])

  tiempo <- data$tiempo
  evento <- data$evento

  pos_pat <- rownames(data[which(data$evento == 1),])
  neg_pat <- rownames(data[which(data$evento == 0),])

  for(x in 1:100) {
    id_pos <- sample(pos_pat, 1)
    id_neg <- sample(neg_pat, 1)
    #test <- as.numeric(c(id_pos, id_neg))
    test <- c(id_pos, id_neg)

    expresion <- as.data.frame(expresion)
    expresion_train <- expresion[-which(rownames(expresion) %in% test),]
    tiempo_train <- data[-which(rownames(expresion) %in% test), "tiempo"]
    evento_train <- data[-which(rownames(expresion) %in% test), "evento"]

    #tiempo_train <- data$tiempo[-test]
    #evento_train <- data$evento[-test]

    expresion_test <- expresion[test,]

    surv_obj <- Surv(tiempo_train, evento_train)
    cox_lasso <- cv.glmnet(as.matrix(expresion_train), surv_obj, family = "cox", alpha = 1)
    best_lam <- cox_lasso$lambda.min
    lasso_best <- suppressWarnings(glmnet(expresion_train, surv_obj, family = "cox", alpha = 1, lambda = best_lam))

    model_matrix <- as.matrix(lasso_best[["beta"]])

    riesgo_train <- predict(lasso_best, newx = as.matrix(expresion_train), type = "link")
    umbral <- median(riesgo_train)

    riesgo_lineal <- predict(lasso_best, newx = as.matrix(expresion_test), type = "link")
    grupo_riesgo <- ifelse(riesgo_lineal > umbral, "High Risk", "Low Risk")

    #grupo_riesgo <- cbind(METADATA_Actualizada[which(METADATA_Actualizada$ID.ENVIRON %in% test), c("ID.ENVIRON", "tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- cbind(metadata[which(rownames(metadata) %in% test), c("subtipo","tiempo", "evento")], grupo_riesgo)
    grupo_riesgo <- as.data.frame(grupo_riesgo)

    grupo_riesgo$predicho <- ifelse(grupo_riesgo$s0 == "High Risk", 1, 0)


    tab <- table(grupo_riesgo$evento, grupo_riesgo$s0)
    stat_iter$tp[x] <- tab[2]
    stat_iter$tn[x] <- tab[3]
    stat_iter$fn[x] <- tab[4]
    stat_iter$fp[x] <- tab[1]
    list_model[[x]] <- grupo_riesgo
    list_genes[[x]] <- model_matrix
  }

  stat_iter[is.na(stat_iter)] <- 0
  list_genes <- do.call(cbind, list_genes)
  colnames(list_genes) <- 1:100
  result <- list(stat_iter, list_model, list_genes)
  return(result)
}

