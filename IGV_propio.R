# ... tu código anterior ...

# Dentro de tu bucle 'for' (asumiendo que 'fusion' es tu iterador actual)
# Reemplaza la llamada a readCoverage por readAlignmentsForAnalysis para obtener las alineaciones completas

# Primero, obtener las regiones de cobertura como ya lo haces
coverageRegion1 <- determineCoverageRegion(exons, fusions[fusion,"gene_id1"], fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], showVicinity[1], showVicinity[2])
coverageRegion2 <- determineCoverageRegion(exons, fusions[fusion,"gene_id2"], fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], showVicinity[3], showVicinity[4])

# Cargar las alineaciones completas para el análisis
alignments1_full <- readAlignmentsForAnalysis(alignmentsFile, fusions[fusion,"contig1"], coverageRegion1)
alignments2_full <- readAlignmentsForAnalysis(alignmentsFile, fusions[fusion,"contig2"], coverageRegion2)

# Ahora puedes calcular la cobertura a partir de estas alineaciones completas
coverage1 <- coverage(alignments1_full)[[as.character(fusions[fusion,"contig1"])]]
coverage2 <- coverage(alignments2_full)[[as.character(fusions[fusion,"contig2"])]]

if (length(alignments1_full) > 0) {
  correct_contig_name1 <- as.character(seqnames(alignments1_full)[1])
  coverage1 <- coverage(alignments1_full)[[correct_contig_name1]]
} else {
  coverage1 <- Rle(0, 0) # Si no hay alineaciones, la cobertura es 0
}

if (length(alignments2_full) > 0) {
  correct_contig_name2 <- as.character(seqnames(alignments2_full)[1])
  coverage2 <- coverage(alignments2_full)[[correct_contig_name2]]
} else {
  coverage2 <- Rle(0, 0)
}

# Ahora, revisa esta parte de tu código, ya que también puede tener problemas si coverage1 es NULL o vacío
# shrink coverage range to chromosome boundaries to avoid subscript out of bounds errors
if (!is.null(coverage1) && length(coverage1) > 0) {
  coverageRegion1 <- IRanges(max(start(coverageRegion1), 1), min(end(coverageRegion1), length(coverage1)))
} else {
  # Si no hay cobertura, ajusta coverageRegion1 para reflejarlo, o maneja el caso
  # (por ejemplo, establece una región vacía o un valor por defecto)
  coverageRegion1 <- IRanges(1,1) # O IRanges() para un rango vacío
  warning(paste("No coverage data found for contig", fusions[fusion,"contig1"]))
}
if (!is.null(coverage2) && length(coverage2) > 0) {
  coverageRegion2 <- IRanges(max(start(coverageRegion2), 1), min(end(coverageRegion2), length(coverage2)))
} else {
  coverageRegion2 <- IRanges(1,1)
  warning(paste("No coverage data found for contig", fusions[fusion,"contig2"]))
}


# ... el resto de tu código de normalización de cobertura y manejo de exons ...

# *************************************************************************
# AHORA, LA PARTE PARA LAS SPLIT READS (como la definimos antes)
# *************************************************************************

# Definir una ventana de proximidad para el breakpoint.
vicinity_for_split_reads <- 20

# Contar split reads para el primer breakpoint
split_reads_count1 <- countSplitReadsAtBreakpoint(
  alignments = alignments1_full,
  breakpoint = fusions[fusion,"breakpoint1"],
  vicinity_bp = vicinity_for_split_reads
)
cat(paste0("Fusión ", fusion, " - Gene1 (", fusions[fusion,"gene_id1"], ") - Breakpoint ", fusions[fusion,"breakpoint1"], ": ",
           split_reads_count1, " lecturas divididas.\n"))

# Contar split reads para el segundo breakpoint
split_reads_count2 <- countSplitReadsAtBreakpoint(
  alignments = alignments2_full,
  breakpoint = fusions[fusion,"breakpoint2"],
  vicinity_bp = vicinity_for_split_reads
)
cat(paste0("Fusión ", fusion, " - Gene2 (", fusions[fusion,"gene_id2"], ") - Breakpoint ", fusions[fusion,"breakpoint2"], ": ",
           split_reads_count2, " lecturas divididas.\n"))

###########################################################################################

# Función para añadir 'chr' si es necesario
addChr <- function(contig_name) {
  if (!grepl("^chr", contig_name, ignore.case = TRUE)) {
    return(paste0("chr", contig_name))
  }
  return(contig_name)
}

library(GenomicAlignments)
library(Rsamtools)
library(GenomicRanges)
library(IRanges) # Asegúrate de que todas estas librerías estén cargadas

# Función para leer alineaciones completas con información de CIGAR
# Esta función reemplazará la forma en que cargas 'alignments' dentro de tu bloque 'if (alignmentsFile != "")'
readAlignmentsForAnalysis <- function(alignmentsFile, contig, region) {
  # Definimos qué campos del BAM queremos leer.
  # 'cigar' es crucial para las split reads.
  # 'qname' (query name) es útil para identificar reads individuales.
  param <- ScanBamParam(
    which = GRanges(contig, region),
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "isize", "mpos", "mrnm"),
    # Opcional: filtrar reads duplicadas o de baja calidad si lo deseas
    flag = scanBamFlag(isDuplicate = FALSE, isNotPassingQualityControls = FALSE)
  )
  
  alignments <- tryCatch(
    {
      readGAlignments(alignmentsFile, param=param)
    },
    error=function(e) {
      # Si falla, intenta añadiendo "chr" al nombre del contig
      param_chr <- ScanBamParam(
        which = GRanges(addChr(contig), region),
        what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "isize", "mpos", "mrnm"),
        flag = scanBamFlag(isDuplicate = FALSE, isNotPassingQualityControls = FALSE)
      )
      readGAlignments(alignmentsFile, param=param_chr)
    }
  )
  return(alignments)
}


# Función para contar split reads cerca de un breakpoint
# Función para contar split reads cerca de un breakpoint
split_reads_count2 <- countSplitReadsAtBreakpoint(
  alignments = alignments2_full,
  breakpoint = fusions[fusion,"breakpoint2"],
  vicinity_bp = vicinity_for_split_reads
)
countSplitReadsAtBreakpoint <- function(alignments, breakpoint, vicinity_bp=10) {
  if (length(alignments) == 0) return(0)
  
  # Filtrar alineaciones que tienen soft-clipping (S en la cadena CIGAR)
  has_soft_clip <- grepl("S", cigar(alignments))
  soft_clipped_reads <- alignments[has_soft_clip]
  
  if (length(soft_clipped_reads) == 0) return(0)
  
  count <- 0
  for (i in seq_along(soft_clipped_reads)) {
    al <- soft_clipped_reads[i]
    
    # *******************************************************************
    # CORRECCIÓN AQUÍ: Asegurarse de que read_cigar sea una cadena de caracteres simple
    read_cigar <- as.character(cigar(al)) 
    # *******************************************************************
    
    read_start <- start(al)
    
    # Obtener las operaciones CIGAR y sus longitudes
    # Esto ya debería funcionar correctamente con una cadena simple
    ops_list <- GenomicAlignments::explodeCigarOps(read_cigar)
    lengths_list <- GenomicAlignments::explodeCigarOpLengths(read_cigar)
    
    # Asegurarnos de que 'ops' y 'lengths' sean vectores atómicos (un solo CIGAR)
    # Como estamos iterando sobre 'soft_clipped_reads' y tomando un elemento a la vez,
    # 'ops_list' y 'lengths_list' serán listas de un solo elemento.
    # Necesitamos extraer el vector de dentro de esa lista.
    ops <- ops_list[[1]]      
    lengths <- lengths_list[[1]]
    
    # Calcular la posición del final alineado de la lectura
    # Excluye S (soft clip), H (hard clip), P (padding)
    aligned_length <- sum(lengths[ops %in% c("M", "I", "D", "N", "EQ", "X")])
    aligned_end <- read_start + aligned_length - 1 # Posición final en la referencia (1-based)
    
    # Caso 1: Soft clip al inicio (S al principio del CIGAR)
    if (ops[1] == "S") {
      if (abs(aligned_end - breakpoint) <= vicinity_bp) {
        count <- count + 1
      }
    # Caso 2: Soft clip al final (S al final del CIGAR)
    } else if (ops[length(ops)] == "S") {
      if (abs(read_start - breakpoint) <= vicinity_bp) {
        count <- count + 1
      }
    }
  }
  return(count)
}

###########################################################################################
###########################################################################################

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

# EXONS ------------------------------------------------------------------------------------------
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


