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

plot_cnv <- function(CNV_calls) {

  data("exons.hg19")
  cnv <- subset(CNV_calls, select = c("start", "end", "name.gene", "type"))
  colnames(cnv)[3] <- "name"
  cnv$col <- ifelse(cnv$type == "duplication", "green", "red")
  bed$col <- "black"

  todo <- rbind(bed[, 2:5], cnv[, c(1, 2, 3, 5)])
  todo <- todo[order(todo$start), ]
  todo$yname <- ifelse(todo$col == "black", ifelse(seq_along(todo$col) %% 2 == 0, 0.25, -0.25), ifelse(todo$col == "green", 0.75, -0.75))
  todo$y <- ifelse(todo$col == "black", 0, ifelse(todo$col == "green", 1, -1))


  #Ordenado segun start:
  gene_order <- unique(todo[order(todo$start), ]$name)
  todo$name <- factor(todo$name, levels = gene_order)
  #todo$col <- factor(todo$col)

  #Delete the rows which names are duplicated but keeping the one that is coloured:
  todo <- ungroup(filter(group_by(todo, name), !(n() > 1 & col == "black")))

  showtext_auto()
  graph <- ggplot(todo, aes(x = start, y = name)) +
    geom_segment(aes(xend = end, yend = name), color = todo$col,  size = 1.2) +
    scale_color_manual(values = c("red", "green", "black")) +
    scale_fill_manual(values = c("red", "green", "black")) +
    labs(x = "Posición", y = "Genes", color = "Color", fill = "Color") +
    theme_classic() +
    theme(axis.text.y = element_markdown(color = todo$col),
          axis.text.y.right = element_blank(),
          axis.text.y.left = element_text(vjust = 1, size= 30),
          axis.title = element_text(size = 50),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(size = 50) )
  graph
  #plot_cnv<- recordPlot()
  return(graph)
}

#' @title Get CNVs annotation.
#' @description Generates a graph that represents the CNVs detected. The user will be able to
#' visualize the distribution, sizes, and positions of the variations within the whole mtDNA
#' and the genes that are affected by them.
#' @param chr dataframe that contains the information about the CNVs detected.
#' @param countThreshold is a number
#' @return lckerm wdkme
#' @export
#'
#x <- all.exons
#chr <- 1
GetCNVsAnnotation <- function(chr, x, countThreshold =10){
  anno <- x@annotations
  selected <- which(anno$chromosome == chr & (x@test + x@reference) *
                      x@expected > countThreshold)

  if (length(selected) == 0) {
    warning("No exon seem to be located in the requested region. It could also be that the read count is too low for these exons? In any case no graph will be plotted.")
    return(NULL)
  }
  anno <- anno[selected, ]
  anno$expected <- x@expected[selected]
  anno$freq <- x@test[selected]/(x@reference[selected] +
                                   x@test[selected])
  anno$middle <- 0.5 * (anno$start + anno$end)
  anno$ratio <- anno$freq/anno$expected
  anno$test <- x@test[selected]
  anno$reference <- x@reference[selected]
  anno$total.counts <- anno$test + anno$reference

  if (length(x@phi) == 1) {
    anno$phi <- x@phi
  } else {
    anno$phi <- x@phi[selected]
  }
  anno <- cbind(anno,plyr::ldply(1:nrow(anno), function(i){
    c(my.min.norm.prop=qbetabinom(p = 0.025, size = anno$total.counts[i],
                                  phi = anno$phi[i], prob = anno$expected[i]),
      my.max.norm.prop=qbetabinom(p = 0.975, size = anno$total.counts[i],
                                  phi = anno$phi[i], prob = anno$expected[i]) )
  }))

  anno$my.min.norm.prop <- anno$my.min.norm/anno$total.counts
  anno$my.max.norm.prop <- anno$my.max.norm/anno$total.counts
  CNVsret <- subset(x@CNV.calls, chromosome == chr)
  attr(CNVsret,"BamHeader") <- attr(x@CNV.calls,"BamHeader")
  return(list(Exons=anno, CNVs=CNVsret))
}


GetAllExons <- function(cnvCalls = all.exons, countThreshold =10){
  #chrs <- unique(cnvCalls$CNVs@CNV.calls$chromosome)
  chrs <- unique(cnvCalls@CNV.calls$chromosome)
  chr.exons <- bplapply(chrs, GetCNVsAnnotation, x=cnvCalls, BPPARAM = bpparam())
  names(chr.exons) <- chrs
  return(invisible(chr.exons))
}

.PlotCNVchromosome <- function(annot = all.exons, thLength=500, geneTable){
  ##assume chrome.bed
  # annot <- chr19
  # thLength <- 500
  # annot <- CHRi

  #bh <- attr(annot$CNVs,"BamHeader")
  #bh <- attr(annot,"BamHeader")

  #data(Centromere)
  load("/media/16TBDisk/Daniela/OMICsdo/data/Centromere.RData")

  Genome <- "hg19"

  #if(!is.null(bh$GenomeDBversion)){
  #  Genome <- ifelse(stringr::str_detect(bh$GenomeDBversion,"GRCh38"),"GRCh38","hg19")
  #}else{
  #  Genome <- ifelse(!any(stringr::str_detect(bh$Code,"hg19")),"GRCh38","hg19")
  #}


  annot <- GetCNVsAnnotation(chr = 1, x = all.exons)
  CNVs <- annot$CNVs
  CNVs$width <- (CNVs$end-CNVs$start)
  chr <- stringr::str_replace_all(annot$Exons$chromosome[1],"chr","Chr")
  CNVs <- subset(CNVs, width > thLength)

  annot <- annot$Exons
  ordm <- order(annot$middle)
  puntos <- data.frame(x=annot$middle[ordm], y=annot$ratio[ordm])
  puntos$colscale <- puntos$y
  puntos$colscale[puntos$colscale>2] <- 2
  # puntos$Type <- NA
  # puntos$Type[puntos$y > annot$my.max.norm.prop] <- "Gain"
  # puntos$Type[puntos$y < annot$my.min.norm.prop] <- "Loss"
  puntos$ls <- annot$my.max.norm.prop[ordm]/annot$expected[ordm]
  puntos$li <- annot$my.min.norm.prop[ordm]/annot$expected[ordm]
  # lsup <- loess(ls~x,puntos, span = 0.3)
  # linf <- loess(li ~x, puntos, span = 0.3)
  # df.l<- data.frame(x=puntos$x,y=lsup$fitted)
  # df.li<- data.frame(x=puntos$x,y=linf$fitted)
  chr <- unique(annot$chromosome)
  if (nrow(subset(puntos,x < Centromere[[Genome]][Centromere[[Genome]]$chr==chr,]$left))>0) {
    lp <- loess(y~x,subset(puntos,x < Centromere[[Genome]][Centromere[[Genome]]$chr==chr,]$left), span = 0.5)
  } else {
    lp<-data.frame(x=puntos$x[1],y=puntos$y[1],fitted=puntos$x[1])
  }


  lq <- loess(y~x,subset(puntos,x > Centromere[[Genome]][Centromere[[Genome]]$chr==chr,]$left), span = 0.5)

  p <- ggplot(puntos, aes(x,y,colour=colscale)) + geom_point(size=0.5) +
    scale_color_stepsn(colours = c("red", "white",  "blue"),breaks = c(  0.9, 1.1 ) )+
    geom_line(aes(x=x,y=y),data=data.frame(x=lp$x,y=lp$fitted),colour="red")+
    geom_line(aes(x=x,y=y),data=data.frame(x=lq$x,y=lq$fitted),colour="red")+
    theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=1, size=0.5)) +
    ylim(-.1,3.5) + geom_line(aes(x=x,y=y),data=data.frame(x=lp$x,y=lp$fitted),colour="red")+
    geom_line(aes(x=x,y=y),data=data.frame(x=lq$x,y=lq$fitted),colour="red")+
    geom_hline(yintercept = 1)

  # p <- ggplot(puntos, aes(x,y,colour=colscale)) + geom_point(size=0.5) +
  #   scale_color_stepsn(colours = c("red", "white",  "blue"),breaks = c(  0.9, 1.1 ) ) +
  #   geom_smooth(method="loess",span=0.1, se=F,col="black", size=0.5,linetype = "dashed") +
  #   # geom_smooth(data = puntos,  method="loess",span=0.1, se=F, aes(x=x,y=ls),col="aquamarine4", size=0.5,linetype = "dashed") +
  # geom_smooth(data = puntos,  method="loess",span=0.1, se=F, aes(x,y=li),col="aquamarine4", size=0.5,linetype = "dashed") +
  # theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=1, size=0.5)) + ylim(-.1,3.5) + geom_hline(yintercept = 1)

  # geom_vline(xintercept = chrom.bed[chr,"mid.point"], linetype = "dashed" ) + xlab("Chr. position") + ylab("log ratio")
  if(nrow(CNVs)>0){
    if(any(CNVs$reads.ratio>1)){
      p <- p + geom_segment(data=subset(CNVs,reads.ratio>1), aes(x = start, y = reads.ratio, xend = end, yend = reads.ratio), col="blue", size=2)
    }
    if(any(CNVs$reads.ratio<1)){
      p <- p + geom_segment(data=subset(CNVs,reads.ratio<1), aes(x = start, y = reads.ratio, xend = end, yend = reads.ratio), col="red", size=2)
    }
  }
  if(missing(geneTable) == FALSE){
    chr <- unique(annot$chromosome)
    genes <- subset(geneTable, chromosome == chr)
    if(nrow(genes)>0){
      exons <- plyr::ldply(genes$Gene,function(x){
        ret <- subset(annot, stringr::str_detect(name, paste0("^",x,"_")))
        if(nrow(ret)!=0){
          ret$gene <- x
          return(ret)
        }

      })
      if(nrow(exons)>0){
        require(dplyr)
        require(ggrepel)
        labels <- exons %>% group_by(gene) %>% summarise(x=min(middle),y=max(ratio))
        p <- p + geom_point(data=exons, mapping=aes(x=middle,y=ratio, shape=gene), colour = "black") +
          geom_text_repel(data=as.data.frame(labels), aes(x,y=3, label=gene, size=0.6), hjust=0, vjust=0,angle=90,colour="black")
      }

    }

  }
  p <- p + labs(title = paste0("Chr ",unique(annot$chromosome)),x="base position",y="Ratio") +
    theme(legend.position = "none",plot.title = element_text(margin = margin(t=40,b=-30)))
  print(p)

  return(invisible(p))
}
