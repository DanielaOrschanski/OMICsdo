#' @title Generation of the plot for CNVs.
#' @description Generates a graph that represents the CNVs detected. The user will be able to
#' visualize the distribution, sizes, and positions of the variations within the whole mtDNA
#' and the genes that are affected by them.
#' @param CNV_calls dataframe that contains the information about the CNVs detected.
#' @return graph
#' @export
#' @import ggplot2
#' @import tidyverse
#' @import ggtext
#' @import dplyr
#' @import showtext
#' @import tidyr

plot_cnv <- function(all.exons, patient_dir, Nspan = 0.8) {

  library(ggplot2)
  library(dplyr)

  x <- all.exons
  anno <- all.exons@annotations
  anno$expected <- x@expected
  anno$freq <- x@test/(x@reference + x@test)
  anno$mitad <- 0.5 * (anno$start + anno$end)
  anno$reads.ratio <- anno$freq/anno$expected
  anno$test <- x@test
  anno$reference <- x@reference
  anno$total.counts <- anno$test + anno$reference
  colnames(anno)[1] <- "Exons"


  expanded_df <- anno

  expanded_df$color <- ifelse(expanded_df$reads.ratio > 1, "duplication", "deletion")
  expanded_df$color <- ifelse(expanded_df$reads.ratio == 1, "no CNV", expanded_df$color)


  # Agrego el maximo anterior a cada valor de start, end y mitad para que queden ordenados por cromosoma ------------------------------
  cromosomas <- 1:22
  expanded_df$globalpos <- NA
  str(expanded_df)
  expanded_df$globalpos[which(expanded_df$chromosome == cromosomas[1])] <- expanded_df$mitad[expanded_df$chromosome == cromosomas[1]]
  expanded_df$globalstart <- NA
  expanded_df$globalstart[which(expanded_df$chromosome == cromosomas[1])] <- expanded_df$start[expanded_df$chromosome == cromosomas[1]]
  expanded_df$globalend <- NA
  expanded_df$globalend[which(expanded_df$chromosome == cromosomas[1])] <- expanded_df$end[expanded_df$chromosome == cromosomas[1]]

  load("/media/16TBDisk/Daniela/OMICsdo/data/Centromere.RData")
  c <- Centromere$hg19
  c$global_left <- NA
  c$global_right <- NA
  c$global_left[which(c$chr == cromosomas[1])] <- c$left[which(c$chr == cromosomas[1])]
  c$global_right[which(c$chr == cromosomas[1])] <- c$right[which(c$chr == cromosomas[1])]

  i=2
  for(i in 2:22) {
    print(cromosomas[i])
    max_anterior <-  max(expanded_df$globalpos[expanded_df$chromosome == cromosomas[i-1]])
    expanded_df$globalpos[which(expanded_df$chromosome == cromosomas[i])] <- expanded_df$mitad[expanded_df$chromosome == cromosomas[i]] + max_anterior
    expanded_df$globalstart[which(expanded_df$chromosome == cromosomas[i])] <- expanded_df$start[expanded_df$chromosome == cromosomas[i]] + max_anterior
    expanded_df$globalend[which(expanded_df$chromosome == cromosomas[i])] <- expanded_df$end[expanded_df$chromosome == cromosomas[i]] + max_anterior

    c$global_left[which(c$chr == cromosomas[i])] <- c$left[which(c$chr == cromosomas[i])]  + max_anterior
    c$global_right[which(c$chr == cromosomas[i])] <- c$right[which(c$chr == cromosomas[i])]  + max_anterior

    i <- i+1
  }

  str(expanded_df$chromosome)


  expanded_df <- expanded_df %>%
    mutate(globalpos = as.numeric(globalpos))

  #Extraigo las posiciones máximas de start y end entre los cromosomas para poner las lineas verticales ---------------------------
  max_positions <- data.frame(
    chromosome = unique(expanded_df$chromosome),
    max_globalpos = sapply(unique(expanded_df$chromosome), function(x) {
      max(expanded_df$globalpos[expanded_df$chromosome == x], na.rm = TRUE)
    }),
    max_globalend = sapply(unique(expanded_df$chromosome), function(x) {
      max(expanded_df$globalend[expanded_df$chromosome == x], na.rm = TRUE)
    }),
    min_globalstart = sapply(unique(expanded_df$chromosome), function(x) {
      min(expanded_df$globalstart[expanded_df$chromosome == x], na.rm = TRUE)
    })
  )


  # Genero una línea que marque el promedio de los ratios por cada brazo de cada cromosoma --------------------------------------------
  medias_p <- data.frame(chr = 1:22, media_p = 0, start = 0, end = 0)
  medias_q <- data.frame(chr = 1:22, media_q = 0, start = 0, end = 0)
  list_loes_p <- list()
  list_loes_q <- list()

  for (i in 1:22) {
    print(i)
    expanded_chr <- expanded_df[which(expanded_df$chromosome == i),]

    if (any(is.na(expanded_chr$reads.ratio))) {
      expanded_chr <- expanded_chr[-which(is.na(expanded_chr$reads.ratio)),]
    }

    centromere_p <- c$right[which(c$chr == i)]

    #chr13 <- bed_axel[which(bed_axel$chromosome == i),]
    #chr13 <- exons.hg19[which(exons.hg19$chromosome == i),]
    #any(chr13$start < centromere_p)

    #Valores medios -----------------------------------------------------------------------------
    mean_p <- mean(expanded_chr$reads.ratio[which(expanded_chr$mitad < centromere_p)])
    mean_q <- mean(expanded_chr$reads.ratio[which(expanded_chr$mitad > centromere_p)])

    medias_p[i,"chr"] <- i
    medias_p[i, "media_p"] <- mean_p
    medias_p[i, "start"] <- min(expanded_chr$globalstart)
    medias_p[i, "end"] <- c$global_left[which(c$chr == i)]

    medias_q[i,"chr"] <- i
    medias_q[i, "media_q"] <- mean_q
    medias_q[i, "start"] <- c$global_left[which(c$chr == i)]
    medias_q[i, "end"] <- max(expanded_chr$globalend)

    #loes ----------------------------------------------------
    # Datos para el brazo p
    p_arm <- expanded_chr[which(expanded_chr$mitad < centromere_p), ]
    if(nrow(p_arm) != 0) {
      loess_p <- loess(reads.ratio ~ globalpos, data = p_arm, span = Nspan)
      loes_p_df <- data.frame(x= loess_p$x, fitted = loess_p$fitted)
    } else {
      loes_p_df <- data.frame()
    }

    # Datos para el brazo q
    q_arm <- expanded_chr[which(expanded_chr$mitad > centromere_p), ]
    loess_q <- loess(reads.ratio ~ globalpos, data = q_arm, span = Nspan)
    loes_q_df <- data.frame(x= loess_q$x, fitted = loess_q$fitted)

    list_loes_p[[i]] <- loes_p_df
    list_loes_q[[i]] <- loes_q_df
  }

  medias_p[is.na(medias_p)] <- 1
  medias_q[is.na(medias_q)] <- 1


  c$chr <- paste("chr", c$chr, sep=" ")

  nas <- expanded_df[which(is.na(expanded_df$reads.ratio)),]
  expanded_df$reads.ratio[which(is.na(expanded_df$reads.ratio))] <- 1
  expanded_df$color[which(expanded_df$reads.ratio == 1)] <- "no CNV"
  expanded_df$colorcnvs <- ifelse(expanded_df$Exons %in% exons_cnvs$Exons, "CNV", "no CNV")
  length(which(expanded_df$colorcnvs == "CNV"))

  g <- ggplot(expanded_df, aes(x = globalpos, y = (reads.ratio), color = color)) +
    geom_point(size= 0.5) +
    geom_segment(data = medias_p, aes(x = start, y = (media_p), xend = end, yend = (media_p)), col="orange", size=1) +
    geom_segment(data = medias_q, aes(x = start, y = (media_q), xend = end, yend = (media_q)), col="orange", size=1) +
    scale_color_manual(values = c("purple", "lightblue", "black")) +
    theme_minimal() +
    labs(x = "Posición dentro del cromosoma (mitad)", y = "Reads Ratio", title = "Distribución de reads.ratio por cromosoma") +
    #geom_vline(data = max_positions, aes(xintercept = min_globalstart), color = "yellow") +
    geom_vline(xintercept = max_positions$min_globalstart[1], color = "black") +
    geom_vline(data = max_positions, aes(xintercept = max_globalend), color = "black") +
    geom_hline(yintercept = 1, color = "black", linetype = "solid", size = 1) +
    scale_x_continuous(breaks = c$global_left, labels = c$chr) +
    geom_vline(data = c, aes(xintercept = global_left), linetype = "dotted", color = "black") +
    scale_y_continuous(breaks = seq(-5, 5, by = 0.5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
  print(g)

  g2 <- g
  for (i in 1:22) {
    g2 <- g2 + geom_line(aes(x = globalpos, y = fitted), data = list_loes_q[[i]], colour = "red", size = 0.3)

    if (nrow(list_loes_p[[i]]) != 0) {
      g2 <- g2 + geom_line(aes(x = globalpos, y = fitted), data = list_loes_p[[i]], colour = "red", size = 0.3)
    }
  }

  print(g2)
  ggsave(filename = sprintf("%s/CNV_graph.png", patient_dir), plot = g2, width = 7, height = 4, dpi = 250)

  return(g2)

}
