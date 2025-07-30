# main loop starts here
#fusion = 2
#alignmentsFile = "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_sorted.bam"

for (fusion in 1:nrow(fusions)) {

  message(paste0("Drawing fusion #", fusion, ": ", fusions[fusion,"gene1"], ":", fusions[fusion,"gene2"]))

  # if showIntergenicVicinity is a number, take it as is
  # if it is a keyword (closestGene/closestProteinCodingGene), determine the range dynamically
  showVicinity <- rep(0, 4)
  if (fusions[fusion,"site1"] == "intergenic") {
    showVicinity[1] <- ifelse(
      is.numeric(showIntergenicVicinity[[1]]),
      showIntergenicVicinity[[1]],
      fusions[fusion,"breakpoint1"] - start(findClosestGene(exons, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], exons$end < fusions[fusion,"breakpoint1"] & exons$type == showIntergenicVicinity[[1]]))
    )
    showVicinity[2] <- ifelse(
      is.numeric(showIntergenicVicinity[[2]]),
      showIntergenicVicinity[[2]],
      end(findClosestGene(exons, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], exons$start > fusions[fusion,"breakpoint1"] & exons$type == showIntergenicVicinity[[2]])) - fusions[fusion,"breakpoint1"]
    )
  }
  if (fusions[fusion,"site2"] == "intergenic") {
    showVicinity[3] <- ifelse(
      is.numeric(showIntergenicVicinity[[3]]),
      showIntergenicVicinity[[3]],
      fusions[fusion,"breakpoint2"] - start(findClosestGene(exons, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], exons$end < fusions[fusion,"breakpoint2"] & exons$type == showIntergenicVicinity[[3]]))
    )
    showVicinity[4] <- ifelse(
      is.numeric(showIntergenicVicinity[[4]]),
      showIntergenicVicinity[[4]],
      end(findClosestGene(exons, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], exons$start > fusions[fusion,"breakpoint2"] & exons$type == showIntergenicVicinity[[4]])) - fusions[fusion,"breakpoint2"]
    )
  }

  # compute coverage from alignments file
  coverage1 <- NULL
  coverage2 <- NULL
  if (alignmentsFile != "") {
    # determine range in which we need to compute the coverage
    #library(IRanges)
    determineCoverageRegion <- function(exons, geneID, contig, breakpoint, showVicinityLeft, showVicinityRight) {
      closestGene <- findClosestGene(exons, contig, breakpoint, exons$geneID == geneID)
      return(IRanges(min(start(closestGene), breakpoint-showVicinityLeft), max(end(closestGene), breakpoint+showVicinityRight)))
    }
    coverageRegion1 <- determineCoverageRegion(exons, fusions[fusion,"gene_id1"], fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], showVicinity[1], showVicinity[2])
    coverageRegion2 <- determineCoverageRegion(exons, fusions[fusion,"gene_id2"], fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], showVicinity[3], showVicinity[4])
    # function which reads alignments from BAM file with & without "chr" prefix
    readCoverage <- function(alignmentsFile, contig, coverageRegion) {
      coverageData <- tryCatch(
        {
          alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(contig, coverageRegion)))
          coverage(alignments)[[contig]]
        },
        error=function(e) {
          alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(addChr(contig), coverageRegion)))
          coverage(alignments)[[addChr(contig)]]
        }
      )
      if (exists("alignments")) rm(alignments)
      return(coverageData)
    }
    # get coverage track
    #library(GenomicAlignments)
    coverage1 <- readCoverage(alignmentsFile, fusions[fusion,"contig1"], coverageRegion1)
    coverage2 <- readCoverage(alignmentsFile, fusions[fusion,"contig2"], coverageRegion2)
    # shrink coverage range to chromosome boundaries to avoid subscript out of bounds errors
    coverageRegion1 <- IRanges(max(start(coverageRegion1), min(start(coverage1))), min(end(coverageRegion1), max(end(coverage1))))
    coverageRegion2 <- IRanges(max(start(coverageRegion2), min(start(coverage2))), min(end(coverageRegion2), max(end(coverage2))))
  }

  # find all exons belonging to the fused genes
  exons1 <- findExons(exons, fusions[fusion,"contig1"], fusions[fusion,"gene_id1"], fusions[fusion,"direction1"], fusions[fusion,"breakpoint1"], coverage1, fusions[fusion,"transcript_id1"], transcriptSelection)
  if (nrow(exons1) == 0) {
    par(mfrow=c(1,1))
    plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
    text(0, 0, paste0("exon coordinates of ", fusions[fusion,"gene1"], " not found in\n", exonsFile))
    warning(paste("exon coordinates of", fusions[fusion,"gene1"], "not found"))
    next
  }
  exons2 <- findExons(exons, fusions[fusion,"contig2"], fusions[fusion,"gene_id2"], fusions[fusion,"direction2"], fusions[fusion,"breakpoint2"], coverage2, fusions[fusion,"transcript_id2"], transcriptSelection)
  if (nrow(exons2) == 0) {
    par(mfrow=c(1,1))
    plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
    text(0, 0, paste0("exon coordinates of ", fusions[fusion,"gene2"], " not found in\n", exonsFile))
    warning(paste("exon coordinates of", fusions[fusion,"gene2"], "not found"))
    next
  }

  # in case of intergenic breakpoints, show the vicinity
  if (sum(showVicinity) > 0) {
    if (fusions[fusion,"site1"] == "intergenic") {
      for (geneID in unique(exons[exons$contig == fusions[fusion,"contig1"] & exons$exonNumber != "intergenic" &
                                  (between(exons$end, fusions[fusion,"breakpoint1"]-showVicinity[1], fusions[fusion,"breakpoint1"]+showVicinity[2]) |
                                   between(exons$start, fusions[fusion,"breakpoint1"]-showVicinity[1], fusions[fusion,"breakpoint1"]+showVicinity[2])),"geneID"]))
        exons1 <- rbind(exons1, findExons(exons, fusions[fusion,"contig1"], geneID, fusions[fusion,"direction1"], fusions[fusion,"breakpoint1"], coverage1, fusions[fusion,"transcript_id1"], transcriptSelection))
      # crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
      exons1 <- exons1[exons1$start >= fusions[fusion,"breakpoint1"]-showVicinity[1] & exons1$end <= fusions[fusion,"breakpoint1"]+showVicinity[2] | exons1$exonNumber == "intergenic",]
    }
    if (fusions[fusion,"site2"] == "intergenic") {
      for (geneID in unique(exons[exons$contig == fusions[fusion,"contig2"] & exons$exonNumber != "intergenic" &
                                  (between(exons$end, fusions[fusion,"breakpoint2"]-showVicinity[3], fusions[fusion,"breakpoint2"]+showVicinity[4]) |
                                   between(exons$start, fusions[fusion,"breakpoint2"]-showVicinity[3], fusions[fusion,"breakpoint2"]+showVicinity[4])),"geneID"]))
        exons2 <- rbind(exons2, findExons(exons, fusions[fusion,"contig2"], geneID, fusions[fusion,"direction2"], fusions[fusion,"breakpoint2"], coverage2, fusions[fusion,"transcript_id2"], transcriptSelection))
      # crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
      exons2 <- exons2[exons2$start >= fusions[fusion,"breakpoint2"]-showVicinity[3] & exons2$end <= fusions[fusion,"breakpoint2"]+showVicinity[4] | exons2$exonNumber == "intergenic",]
    }
  }

  # normalize coverage
  if (alignmentsFile != "") {
    coverageNormalization <- function(coverage, coverageRegion, exons) {
      max(1, ifelse(
        squishIntrons, # => ignore intronic coverage
        max(as.numeric(coverage[IRanges(sapply(exons$start,max,min(start(coverage))),sapply(exons$end,min,max(end(coverage))))])),
        round(quantile(coverage[coverageRegion], 0.9999)) # ignore coverage spikes from read-attracting regions
      ))
    }
    coverageNormalization1 <- ifelse(head(coverageRange,1) == 0, coverageNormalization(coverage1, coverageRegion1, exons1), head(coverageRange,1))
    coverageNormalization2 <- ifelse(tail(coverageRange,1) == 0, coverageNormalization(coverage2, coverageRegion2, exons2), tail(coverageRange,1))
    if (length(coverageRange) == 1 && coverageRange == 0) { # harmonize scales of gene1 and gene2
      coverageNormalization1 <- max(coverageNormalization1, coverageNormalization2)
      coverageNormalization2 <- max(coverageNormalization1, coverageNormalization2)
    }
    coverage1 <- coverage1/coverageNormalization1
    coverage2 <- coverage2/coverageNormalization2
    coverage1[coverage1 > 1] <- 1
    coverage2[coverage2 > 1] <- 1
  }

  # sort coding exons last, such that they are drawn over the border of non-coding exons
  exons1 <- exons1[order(exons1$start, -rank(exons1$type)),]
  exons2 <- exons2[order(exons2$start, -rank(exons2$type)),]

  # insert dummy exons, if breakpoints are outside the gene (e.g., in UTRs)
  # this avoids plotting artifacts
  breakpoint1 <- fusions[fusion,"breakpoint1"]
  breakpoint2 <- fusions[fusion,"breakpoint2"]
  if (breakpoint1 < min(exons1$start)) {
    exons1 <- rbind(c(exons1[1,"contig"], "dummy", max(1,breakpoint1-1000), max(1,breakpoint1-1000), exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""), exons1)
  } else if (breakpoint1 > max(exons1$end)) {
    exons1 <- rbind(exons1, c(exons1[1,"contig"], "dummy", breakpoint1+1000, breakpoint1+1000, exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""))
  }
  if (breakpoint2 < min(exons2$start)) {
    exons2 <- rbind(c(exons2[1,"contig"], "dummy", max(1,breakpoint2-1000), max(1,breakpoint2-1000), exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""), exons2)
  } else if (breakpoint2 > max(exons2$end)) {
    exons2 <- rbind(exons2, c(exons2[1,"contig"], "dummy", breakpoint2+1000, breakpoint2+1000, exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""))
  }
  exons1$start <- as.integer(exons1$start)
  exons1$end <- as.integer(exons1$end)
  exons2$start <- as.integer(exons2$start)
  exons2$end <- as.integer(exons2$end)

  exons1$left <- exons1$start
  exons1$right <- exons1$end
  exons2$left <- exons2$start
  exons2$right <- exons2$end

  squishedIntronSize <- 200
  if (squishIntrons) {
    # hide introns in gene1
    cumulativeIntronLength <- 0
    previousExonEnd <- -squishedIntronSize
    for (exon in 1:nrow(exons1)) {
      if (breakpoint1 > previousExonEnd+1 && breakpoint1 < exons1[exon,"left"])
        breakpoint1 <- (breakpoint1-previousExonEnd) / (exons1[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
      if (exons1[exon,"left"] > previousExonEnd) {
        cumulativeIntronLength <- cumulativeIntronLength + exons1[exon,"left"] - previousExonEnd - squishedIntronSize
        previousExonEnd <- exons1[exon,"right"]
      }
      if (breakpoint1 >= exons1[exon,"left"] && breakpoint1 <= exons1[exon,"right"]+1)
        breakpoint1 <- breakpoint1 - cumulativeIntronLength
      exons1[exon,"left"] <- exons1[exon,"left"] - cumulativeIntronLength
      exons1[exon,"right"] <- exons1[exon,"right"] - cumulativeIntronLength
    }

    # hide introns in gene2
    cumulativeIntronLength <- 0
    previousExonEnd <- -squishedIntronSize
    for (exon in 1:nrow(exons2)) {
      if (breakpoint2 > previousExonEnd+1 && breakpoint2 < exons2[exon,"left"])
        breakpoint2 <- (breakpoint2-previousExonEnd) / (exons2[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
      if (exons2[exon,"left"] > previousExonEnd) {
        cumulativeIntronLength <- cumulativeIntronLength + exons2[exon,"left"] - previousExonEnd - squishedIntronSize
        previousExonEnd <- exons2[exon,"right"]
      }
      if (breakpoint2 >= exons2[exon,"left"] && breakpoint2 <= exons2[exon,"right"]+1)
        breakpoint2 <- breakpoint2 - cumulativeIntronLength
      exons2[exon,"left"] <- exons2[exon,"left"] - cumulativeIntronLength
      exons2[exon,"right"] <- exons2[exon,"right"] - cumulativeIntronLength
    }
  } else { # don't squish introns
    # shift exon coordinates to align the gene to the left border of the plot
    exons1$right <- exons1$right - min(exons1$left)
    breakpoint1 <- breakpoint1 - min(exons1$left)
    exons1$left <- exons1$left - min(exons1$left)
    exons2$right <- exons2$right - min(exons2$left)
    breakpoint2 <- breakpoint2 - min(exons2$left)
    exons2$left <- exons2$left - min(exons2$left)
  }

  # scale exon sizes to fit on page
  scalingFactor <- max(exons1$right) + max(exons2$right)
  if (fixedScale > 0) {
    if (fixedScale >= scalingFactor) {
      scalingFactor <- fixedScale
    } else {
      warning(paste("fallback to automatic scaling, because value for --fixedScale is too small to fit transcripts on canvas (increase it to", scalingFactor, "to avoid this)"))
    }
  }
  exons1$left <- exons1$left / scalingFactor
  exons1$right <- exons1$right / scalingFactor
  exons2$left <- exons2$left / scalingFactor
  exons2$right <- exons2$right / scalingFactor
  breakpoint1 <- breakpoint1 / scalingFactor
  breakpoint2 <- breakpoint2 / scalingFactor

  # shift gene2 to the right border of the page
  gene2Offset <- 1 + 0.05 - max(exons2$right)

  # center fusion horizontally
  fusionOffset1 <- (max(exons1$right)+gene2Offset)/2 - ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)
  fusionOffset2 <- fusionOffset1 + ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)

  ########################################
  # ARRANCA GRAFICO ########################################
  ####################################

  # layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
  layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
  par(mar=c(0, 0, 0, 0))
  plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.0, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

  # vertical coordinates of layers
  ySampleName <- 1.04
  #yIdeograms <- ifelse(alignmentsFile != "", 0.94, 0.84)
  #yBreakpointLabels <- ifelse(alignmentsFile != "", 0.86, 0.76)
  yBreakpointLabels <- ifelse(alignmentsFile != "", 0.76, 0.76)

  #yCoverage <- 0.72
  yCoverage <- 0.4
  #yExons <- 0.67
  yExons <- 0.3
  #yGeneNames <- 0.58
  yGeneNames <- 0.88
  yFusion <- 0.1
  yTranscript <- 0.45
  yScale <- 0.407
  yTrajectoryBreakpointLabels <- yBreakpointLabels - 0.035
  yTrajectoryExonTop <- yExons + 0.03
  yTrajectoryExonBottom <- yExons - 0.055
  yTrajectoryFusion <- yFusion + 0.03

  # print sample name (title of page)
  text(0.5, ySampleName, sampleName, font=2, cex=fontSize*1.5, adj=c(0.5,0))


  # draw gene & transcript names
  #fusion = 2
  if (fusions[fusion,"gene1"] != ".")
    #text(max(exons1$right)/2, yGeneNames, fusions[fusion,"gene1"], font=2, cex=fontSize, adj=c(0.5,0))
    text(0.25, yGeneNames, fusions[fusion,"gene1"], font=2, cex=fontSize, adj=c(0.5,0))
  if (fusions[fusion,"site1"] != "intergenic")
    #text(max(exons1$right)/2, yGeneNames-0.01, head(exons1$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))
    text(0.25, yGeneNames-0.02, sprintf("transcript id = %s", head(exons1$transcript,1)), cex=0.9*fontSize, adj=c(0.5,1))
  if (fusions[fusion,"gene2"] != ".")
    #text(gene2Offset+max(exons2$right)/2, yGeneNames, fusions[fusion,"gene2"], font=2, cex=fontSize, adj=c(0.5,0))
    text(0.75, yGeneNames, fusions[fusion,"gene2"], font=2, cex=fontSize, adj=c(0.5,0))
  if (fusions[fusion,"site2"] != "intergenic")
    #text(gene2Offset+max(exons2$right)/2, yGeneNames-0.01, head(exons2$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))
    text(0.75, yGeneNames-0.02, sprintf("transcript id = %s", head(exons2$transcript,1)), cex=0.9*fontSize, adj=c(0.5,1))

  # if multiple genes in the vicinity are shown, label them
  if (fusions[fusion,"site1"] == "intergenic")
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene & exons1$type != "dummy",]
      if (any(exonsOfGene$type == "exon"))
        text(mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yGeneNames-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
    }
  if (fusions[fusion,"site2"] == "intergenic")
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene & exons2$type != "dummy",]
      if (any(exonsOfGene$type == "exon"))
        text(gene2Offset+mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yGeneNames-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
    }

  # label breakpoints
  #text(breakpoint1+0.01, yBreakpointLabels-0.03, paste0("bp1 = ", fusions[fusion,"display_contig1"], ":", fusions[fusion,"breakpoint1"]), adj=c(1,0), cex=fontSize)
  text(breakpoint1+0.05, yBreakpointLabels, paste0("bp1 = ", fusions[fusion,"display_contig1"], ":", fusions[fusion,"breakpoint1"]), adj=c(1,0), cex=fontSize)
  #text(gene2Offset+breakpoint2-0.01, yBreakpointLabels-0.03, paste0("breakpoint2\n", fusions[fusion,"display_contig2"], ":", fusions[fusion,"breakpoint2"]), adj=c(0,0), cex=fontSize)
  text(gene2Offset+breakpoint2-0.05, yBreakpointLabels - 0.05, paste0("bp2 = ", fusions[fusion,"display_contig2"], ":", fusions[fusion,"breakpoint2"]), adj=c(0,0), cex=fontSize)

  # draw coverage axis
  if (alignmentsFile != "") {
    # left axis (gene1)
    #yCoverage = 0.4
    #lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
    lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.3, yCoverage+0.3))
    text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.9*fontSize)
    text(-0.025, yCoverage+0.3, coverageNormalization1, adj=c(1,0.5), cex=0.9*fontSize)
    text(-0.05, yCoverage+0.2, "Coverage", srt=90, cex=0.9*fontSize, adj=c(1,0.5))

    # right axis (gene2)
    if (length(coverageRange) == 2) { # separate axes for gene1 and gene2
      rightCoverageAxisX <- gene2Offset+max(exons2$right)
      lines(c(rightCoverageAxisX+0.02, rightCoverageAxisX+0.01, rightCoverageAxisX+0.01, rightCoverageAxisX+0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
      text(rightCoverageAxisX+0.025, yCoverage, "0", adj=c(0,0.5), cex=0.9*fontSize)
      text(rightCoverageAxisX+0.025, yCoverage+0.1, coverageNormalization2, adj=c(0,0.5), cex=0.9*fontSize)
      text(rightCoverageAxisX+0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*fontSize, adj=c(1,0.5))
    }

    # plot coverage 1
    rect(min(exons1$left), yCoverage, max(exons1$right), yCoverage+0.3, col="#eeeeee", border=NA)
    if (squishIntrons) {
      for (exon in 1:nrow(exons1))
        if (exons1[exon,"type"] != "CDS") # don't draw coverage twice for coding regions
          drawCoverage(exons1[exon,"left"], exons1[exon,"right"], yCoverage, coverage1, exons1[exon,"start"], exons1[exon,"end"], color1)
    } else {
      drawCoverage(min(exons1$left), max(exons1$right), yCoverage, coverage1, min(exons1$start), max(exons1$end), color1)
    }

    lines(c(breakpoint1, breakpoint1), c(yCoverage, yCoverage+0.3), col = "purple", lty = 2, lwd = 1)

    # plot coverage 2
    rect(gene2Offset+min(exons2$left), yCoverage, gene2Offset+max(exons2$right), yCoverage+0.3, col="#eeeeee", border=NA)
    if (squishIntrons) {
      for (exon in 1:nrow(exons2))
        if (exons2[exon,"type"] != "CDS") # don't draw coverage twice for coding regions
          drawCoverage(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yCoverage, coverage2, exons2[exon,"start"], exons2[exon,"end"], color2)
    } else {
      drawCoverage(gene2Offset+min(exons2$left), gene2Offset+max(exons2$right), yCoverage, coverage2, min(exons2$start), max(exons2$end), color2)
    }
    lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2), c(yCoverage, yCoverage+0.3), col = "darkorange", lty = 2, lwd = 1)

  }

  #Colorear los exones que van a ir a la fusion: -------------------
    #GEN1:
    #1. Definí los extremos según la dirección
  if (fusions[fusion, "direction1"] == "upstream") {
    rect_x_start <- breakpoint1
    rect_x_end <- max(exons1$right)
  } else if (fusions[fusion, "direction1"] == "downstream") {
    rect_x_start <- min(exons1$left)
    rect_x_end <- breakpoint1
  }
    #2. Dibujá el área sombreada
  rect(rect_x_start, yCoverage, rect_x_end, yCoverage+0.3,
       col = adjustcolor("#DCD0FF", alpha.f = 0.5), border = NA)

  #GEN2:
  #1. Definí los extremos según la dirección
  if (fusions[fusion, "direction2"] == "upstream") {
    rect_x_start <- breakpoint2
    rect_x_end <- max(exons2$right)
  } else if (fusions[fusion, "direction2"] == "downstream") {
    rect_x_start <- min(exons2$left)
    rect_x_end <- breakpoint2
  }
  #2. Dibujá el área sombreada
  rect(gene2Offset + rect_x_start, yCoverage, gene2Offset + rect_x_end, yCoverage+0.3,
       col = adjustcolor("#FFDBBB", alpha.f = 0.5), border = NA)


  #-----------------------------------------------------------------------

  # PLOT FUSION ################################

  # strand 1
  #lines(c(min(exons1$left), max(exons1$right)), c(yExons, yExons), col= color1)

  #for (gene in unique(exons1$geneName))
  #  drawStrand(min(exons1[exons1$geneName == gene,"left"]), max(exons1[exons1$geneName == gene,"right"]), yExons, darkColor1, head(exons1[exons1$geneName == gene,"strand"],1))
  #for (exon in 1:nrow(exons1))
  #  drawExon(exons1[exon,"left"], exons1[exon,"right"], yExons, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])

  # plot gene 2
  #lines(c(gene2Offset, gene2Offset+max(exons2$right)), c(yExons, yExons), col=darkColor2)
  #for (gene in unique(exons2$geneName))
  #  drawStrand(gene2Offset+min(exons2[exons2$geneName == gene,"left"]), gene2Offset+max(exons2[exons2$geneName == gene,"right"]), yExons, darkColor2, head(exons2[exons2$geneName == gene,"strand"],1))
  #for (exon in 1:nrow(exons2))
  #  drawExon(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yExons, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])

  yFusion = 0.3
  # plot gene1 of fusion
  if (fusions[fusion,"direction1"] == "downstream") {
    # plot strands
    lines(c(fusionOffset1, fusionOffset1+breakpoint1), c(yFusion, yFusion), col=darkColor1)
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene,]
      if (min(exonsOfGene$start) <= fusions[fusion,"breakpoint1"])
        drawStrand(fusionOffset1+min(exonsOfGene$left), fusionOffset1+min(breakpoint1, max(exonsOfGene$right)), yFusion, col=darkColor1, exonsOfGene$strand[1])
    }
    # plot exons
    for (exon in 1:nrow(exons1))
      if (exons1[exon,"start"] <= fusions[fusion,"breakpoint1"])
        drawExon(fusionOffset1+exons1[exon,"left"], fusionOffset1+min(breakpoint1, exons1[exon,"right"]), yFusion, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
    # plot trajectories
    #lines(c(0, 0, fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
    #lines(c(breakpoint1, breakpoint1, fusionOffset1+breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
  } else if (fusions[fusion,"direction1"] == "upstream") {
    # plot strands
    lines(c(fusionOffset1, fusionOffset2), c(yFusion, yFusion), col=darkColor1)
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene,]
      if (max(exonsOfGene$end+1) >= fusions[fusion,"breakpoint1"])
        drawStrand(fusionOffset2-max(exonsOfGene$right)+breakpoint1, min(fusionOffset2, fusionOffset2-min(exonsOfGene$left)+breakpoint1), yFusion, col=darkColor1, chartr("+-", "-+", exonsOfGene$strand[1]))
    }
    # plot exons
    for (exon in 1:nrow(exons1))
      if (exons1[exon,"end"]+1 >= fusions[fusion,"breakpoint1"])
        drawExon(fusionOffset1+max(exons1$right)-exons1[exon,"right"], min(fusionOffset2, fusionOffset1+max(exons1$right)-exons1[exon,"left"]), yFusion, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
    # plot trajectories
    #lines(c(max(exons1$right), max(exons1$right), fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
    lines(c(breakpoint1, breakpoint1, fusionOffset1+max(exons1$right)-breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
  }

  # plot gene2 of fusion
  if (fusions[fusion,"direction2"] == "downstream") {
    # plot strands
    lines(c(fusionOffset2, fusionOffset2+breakpoint2), c(yFusion, yFusion), col=darkColor2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene,]
      if (min(exonsOfGene$start) <= fusions[fusion,"breakpoint2"])
        drawStrand(max(fusionOffset2, fusionOffset2+breakpoint2-max(exonsOfGene$right)), fusionOffset2+breakpoint2-min(exonsOfGene$left), yFusion, col=darkColor2, chartr("+-", "-+", exonsOfGene$strand[1]))
    }
    # plot exons
    for (exon in 1:nrow(exons2))
      if (exons2[exon,"start"] <= fusions[fusion,"breakpoint2"])
        drawExon(max(fusionOffset2, fusionOffset2+breakpoint2-exons2[exon,"right"]), fusionOffset2+breakpoint2-exons2[exon,"left"], yFusion, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
    # plot trajectories
    #lines(c(gene2Offset, gene2Offset, fusionOffset2+breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
    #lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
  } else if (fusions[fusion,"direction2"] == "upstream") {
    # plot strands
    lines(c(fusionOffset2, fusionOffset2+max(exons2$right)-breakpoint2), c(yFusion, yFusion), col=darkColor2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene,]
      if (max(exonsOfGene$end+1) >= fusions[fusion,"breakpoint2"])
        drawStrand(max(fusionOffset2, fusionOffset2+min(exonsOfGene$left)-breakpoint2), fusionOffset2+max(exonsOfGene$right)-breakpoint2, yFusion, col=darkColor2, exonsOfGene$strand[1])
    }
    # plot exons
    for (exon in 1:nrow(exons2))
      if (exons2[exon,"end"]+1 >= fusions[fusion,"breakpoint2"])
        drawExon(max(fusionOffset2, fusionOffset2+exons2[exon,"left"]-breakpoint2), fusionOffset2+exons2[exon,"right"]-breakpoint2, yFusion, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
    # plot trajectories
    lines(c(gene2Offset+max(exons2$right), gene2Offset+max(exons2$right), fusionOffset2+max(exons2$right)-breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
    #lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
  }


  # draw circos plot
  if (is.null(cytobands) || !("circlize" %in% names(sessionInfo()$otherPkgs)) || !("GenomicRanges" %in% names(sessionInfo()$otherPkgs))) {
    plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
    plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  } else {
    par(mar=c(2,4,0,0), xpd=NA)
    drawCircos(fusion, fusions, cytobands, minConfidenceForCircosPlot, circosColors)
    par(mar=c(0,0,0,0), xpd=F)
  }

  # draw protein domains
  plot(0, 0, type="l", xlim=c(-0.1, 1.1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  par(xpd=NA)
  if (!is.null(proteinDomains) && "GenomicRanges" %in% names(sessionInfo()$otherPkgs))
    drawProteinDomains(fusions[fusion,], exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors)
  par(xpd=F)

  # print statistics about supporting alignments
  #plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  text(0, 0.575, "Supporting Read Count", font=2, adj=c(0,0), cex=fontSize)
  if ("split_reads" %in% colnames(fusions)) { # STAR-Fusion reports split reads from both breakpoints combined
    text(0, 0.525, paste0("Split reads = ", fusions[fusion,"split_reads"], "\n", "Discordant mates = ", fusions[fusion,"discordant_mates"]), adj=c(0,1), cex=fontSize)
  } else { # Arriba reports split reads separately for the two breakpoints
    text(
      0, 0.525,
      paste0(
        "Split reads - bp1 = ", fusions[fusion,"split_reads1"], "\n",
        "Split reads - bp2 = ", fusions[fusion,"split_reads2"], "\n",
        "Discordant mates = ", fusions[fusion,"discordant_mates"]
      ),
      adj=c(0,1), cex=fontSize
    )
  }

}

devNull <- dev.off()
message("Done")
