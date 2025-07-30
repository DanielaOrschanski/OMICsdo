#!/usr/bin/env Rscript

# print warnings as they happen instead of collecting them for after a loop ends
options(warn=1)

# define valid parameters
parameters <- list(
  fusions=list("fusionsFile", "file", "fusions.tsv", T),
  annotation=list("exonsFile", "file", "annotation.gtf", T),
  output=list("outputFile", "string", "output.pdf", T),
  alignments=list("alignmentsFile", "file", "Aligned.sortedByCoord.out.bam"),
  cytobands=list("cytobandsFile", "file", "cytobands.tsv"),
  minConfidenceForCircosPlot=list("minConfidenceForCircosPlot", "string", "medium"),
  proteinDomains=list("proteinDomainsFile", "file", "protein_domains.gff3"),
  sampleName=list("sampleName", "string", ""),
  squishIntrons=list("squishIntrons", "bool", T),
  printExonLabels=list("printExonLabels", "bool", T),
  render3dEffect=list("render3dEffect", "bool", T),
  pdfWidth=list("pdfWidth", "numeric", 11.692),
  pdfHeight=list("pdfHeight", "numeric", 8.267),
  color1=list("color1", "string", "steelblue"),
  color2=list("color2", "string", "indianred"),
  mergeDomainsOverlappingBy=list("mergeDomainsOverlappingBy", "numeric", 0.9),
  optimizeDomainColors=list("optimizeDomainColors", "bool", F),
  fontSize=list("fontSize", "numeric", 1),
  fontFamily=list("fontFamily", "string", "Helvetica"),
  showIntergenicVicinity=list("showIntergenicVicinity", "string", "0"),
  transcriptSelection=list("transcriptSelection", "string", "provided"),
  fixedScale=list("fixedScale", "numeric", 0),
  coverageRange=list("coverageRange", "string", "0")
)


# print help if necessary
args <- commandArgs(trailingOnly=T)

# make sure mandatory arguments are present
for (parameter in names(parameters)){
  if (length(parameters[[parameter]]) > 3 && parameters[[parameter]][[4]])
    if (!any(grepl(paste0("^--", parameter, "="), args), perl=T))
      stop(paste0("Missing mandatory argument: --", parameter))
  # set default values
  assign(parameters[[parameter]][[1]], ifelse(parameters[[parameter]][[2]] == "file", "", parameters[[parameter]][[3]]))

}


# parse command-line parameters
for (arg in args) {
  argName <- sub("=.*", "", sub("^--", "", arg, perl=T), perl=T)
  argValue <- sub("^[^=]*=", "", arg, perl=T)
  if (!(argName %in% names(parameters)) || !grepl("^--", arg, perl=T))
    stop(paste("Unknown parameter:", arg))
  if (parameters[[argName]][[2]] == "bool") {
    if (argValue %in% c("TRUE", "T", "FALSE", "F")) {
      assign(parameters[[argName]][[1]], as.logical(argValue))
    } else {
      stop(paste0("Invalid argument to --", argName))
    }
  } else if (parameters[[argName]][[2]] == "string") {
    assign(parameters[[argName]][[1]], argValue)
  } else if (parameters[[argName]][[2]] == "numeric") {
    if (is.na(suppressWarnings(as.numeric(argValue))))
      stop(paste0("Invalid argument to --", argName))
    assign(parameters[[argName]][[1]], as.numeric(argValue))
  } else if (parameters[[argName]][[2]] == "file") {
    if (file.access(argValue) == -1)
      stop(paste("Cannot read file:", argValue))
    assign(parameters[[argName]][[1]], argValue)
  }
}


library(readxl)
fusions_report <- read_excel("/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_FusionReport.xlsx")

alignmentsFile = "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_sorted.bam"
fusionsFile = "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_fusions.tsv"
cytobandsFile = "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv"
proteinDomainsFile = "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"
sampleName = "SRR6888826"
exonsFile = "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.110.gtf"
color1="purple"
color2="orange"


# validate values of parameters
if (cytobandsFile == "")
  warning("Missing parameter '--cytobands'. No ideograms and circos plots will be drawn.")
if (!(minConfidenceForCircosPlot %in% c("none", "low", "medium", "high")))
  stop("Invalid argument to --minConfidenceForCircosPlot")
showIntergenicVicinity <- as.list(unlist(strsplit(showIntergenicVicinity, ",", fixed=T)))
if (!(length(showIntergenicVicinity) %in% c(1,4)))
  stop(paste0("Invalid argument to --showIntergenicVicinity"))
showIntergenicVicinity <- lapply(showIntergenicVicinity, function(x) {
  if (x == "closestGene") {
    return("exon")
  } else if (x == "closestProteinCodingGene") {
    return("CDS")
  } else if (is.na(suppressWarnings(as.numeric(x))) || as.numeric(x) < 0) {
    stop(paste0("Invalid argument to --showIntergenicVicinity"))
  } else {
    return(as.numeric(x))
  }
})
if (length(showIntergenicVicinity) == 1)
  showIntergenicVicinity <- rep(showIntergenicVicinity, 4)
if (squishIntrons)
  if (any(!is.numeric(unlist(showIntergenicVicinity))) || any(showIntergenicVicinity > 0))
    stop("--squishIntrons must be disabled, when --showIntergenicVicinity is > 0")
if (!(transcriptSelection %in% c("coverage", "provided", "canonical")))
  stop("Invalid argument to --transcriptSelection")
if (fixedScale < 0)
  stop("Invalid argument to --fixedScale")
if (!(fontFamily %in% names(pdfFonts())))
  stop(paste0("Unknown font: ", fontFamily, ". Available fonts: ", paste(names(pdfFonts()), collapse=", ")))
coverageRange <- suppressWarnings(as.numeric(unlist(strsplit(coverageRange, ",", fixed=T))))
if (!(length(coverageRange) %in% 1:2) || any(is.na(coverageRange)) || any(coverageRange < 0))
  stop("Invalid argument to --coverageRange")

# check if required packages are installed
if (!suppressPackageStartupMessages(require(GenomicRanges)))
  warning("Package 'GenomicRanges' is not installed. No protein domains and circos plots will be drawn.")
if (!suppressPackageStartupMessages(require(circlize)))
  warning("Package 'circlize' is not installed. No circos plots will be drawn.")
if (alignmentsFile != "")
  if (!suppressPackageStartupMessages(require(GenomicAlignments)))
    stop("Package 'GenomicAlignments' must be installed when '--alignments' is used")

# define colors
changeColorBrightness <- function(color, delta) {
  rgb(
    min(255,max(0,col2rgb(color)["red",]+delta)),
    min(255,max(0,col2rgb(color)["green",]+delta)),
    min(255,max(0,col2rgb(color)["blue",]+delta)),
    maxColorValue=255
  )
}
getDarkColor <- function(color) { changeColorBrightness(color, -100) }
getBrightColor <- function(color) { changeColorBrightness(color, +190) }
darkColor1 <- getDarkColor(color1)
darkColor2 <- getDarkColor(color2)
circosColors <- c(translocation="#000000", duplication="#00bb00", deletion="#ff0000", inversion="#0000ff")

# convenience functions to add/remove "chr" prefix
addChr <- function(contig) {
  ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
  sub("^chr", "", sub("^chrM", "MT", contig, perl=T), perl=T)
}

# convenience function to check if a value is between two others
between <- function(value, start, end) {
  value >= start & value <= end
}

# read cytoband annotation ------------------------------------
#cytobandsFile = "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv"
cytobands <- NULL

if (cytobandsFile != "") {
  cytobands <- read.table(cytobandsFile, header=T, colClasses=c("character", "numeric", "numeric", "character", "character"))
  cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]
}

# read exon annotation ------------------------------
#exonsFile <- "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.110.gtf"

message("Loading annotation")
exons <- scan(exonsFile, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="#", quote='"', multi.line=F)
attr(exons, "row.names") <- .set_row_names(length(exons[[1]]))
class(exons) <- "data.frame"
exons <- exons[exons$type %in% c("exon","CDS"),c("contig","type","start","end","strand","attributes")]
exons$contig <- removeChr(exons$contig)
parseGtfAttribute <- function(attribute, gtf) {
  parsed <- sub(paste0(".*", attribute, "[ =]([^;]+).*"), "\\1", gtf$attributes, perl=T)
  failedToParse <- parsed == gtf$attributes
  if (any(failedToParse)) {
    warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " record(s)."))
    parsed <- ifelse(failedToParse, "", parsed)
  }
  return(parsed)
}
exons$geneID <- parseGtfAttribute("gene_id", exons)
exons$geneName <- parseGtfAttribute("gene_name", exons)
exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
exons$transcript <- parseGtfAttribute("transcript_id", exons)
exons$exonNumber <- ifelse(rep(printExonLabels, nrow(exons)), parseGtfAttribute("exon_number", exons), "")

# read protein domain annotation
#proteinDomainsFile <- "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"
proteinDomains <- NULL
if (proteinDomainsFile != "") {
  message("Loading protein domains")
  proteinDomains <- scan(proteinDomainsFile, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="", quote="", multi.line=F)
  attr(proteinDomains, "row.names") <- .set_row_names(length(proteinDomains[[1]]))
  class(proteinDomains) <- "data.frame"
  proteinDomains$color <- parseGtfAttribute("color", proteinDomains)
  proteinDomains$proteinDomainName <- sapply(parseGtfAttribute("Name", proteinDomains), URLdecode)
  proteinDomains$proteinDomainID <- parseGtfAttribute("protein_domain_id", proteinDomains)
}

drawCurlyBrace <- function(left, right, top, bottom, tip) {
  smoothness <- 20
  x <- cumsum(exp(-seq(-2.5, 2.5, len=smoothness)^2))
  x <- x/max(x)
  y <- seq(top, bottom, len=smoothness)
  lines(left+(tip-left)+x*(left-tip), y)
  lines(tip+x*(right-tip), y)
}

####################################################################
##### ESPECIFICO DEL ID ########################################
###########################################################

# read fusions -------------------
#fusionsFile <- "/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_fusions.tsv"
#sampleName <- "SRR6888826"
fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")

if (colnames(fusions)[1] == "X.gene1") { # Arriba output
  colnames(fusions)[colnames(fusions) %in% c("X.gene1", "strand1.gene.fusion.", "strand2.gene.fusion.")] <- c("gene1", "strand1", "strand2")
  fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl=T)
  fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl=T)
  fusions$contig1 <- removeChr(fusions$display_contig1)
  fusions$contig2 <- removeChr(fusions$display_contig2)
  fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl=T))
  fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl=T))
  fusions$split_reads1 <- fusions$split_reads1
  fusions$split_reads2 <- fusions$split_reads2
  fusions$type <- sub(".*(translocation|duplication|deletion|inversion).*", "\\1", fusions$type, perl=T)
  fusions$fusion_transcript <- gsub("[()^$]", "", fusions$fusion_transcript)
} else if (colnames(fusions)[1] == "X.FusionName") { # STAR-Fusion
  fusions$gene1 <- sub("\\^.*", "", fusions$LeftGene, perl=T)
  fusions$gene2 <- sub("\\^.*", "", fusions$RightGene, perl=T)
  fusions$strand1 <- sub(".*:(.)$", "\\1/\\1", fusions$LeftBreakpoint, perl=T)
  fusions$strand2 <- sub(".*:(.)$", "\\1/\\1", fusions$RightBreakpoint, perl=T)
  fusions$display_contig1 <- sub(":[^:]*:[^:]*$", "", fusions$LeftBreakpoint, perl=T)
  fusions$display_contig2 <- sub(":[^:]*:[^:]*$", "", fusions$RightBreakpoint, perl=T)
  fusions$contig1 <- removeChr(fusions$display_contig1)
  fusions$contig2 <- removeChr(fusions$display_contig2)
  fusions$breakpoint1 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$LeftBreakpoint, perl=T))
  fusions$breakpoint2 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$RightBreakpoint, perl=T))
  fusions$direction1 <- ifelse(grepl(":\\+$", fusions$LeftBreakpoint, perl=T), "downstream", "upstream")
  fusions$direction2 <- ifelse(grepl(":\\+$", fusions$RightBreakpoint, perl=T), "upstream", "downstream")
  fusions$gene_id1 <- sub(".*\\^", "", fusions$LeftGene, perl=T)
  fusions$gene_id2 <- sub(".*\\^", "", fusions$RightGene, perl=T)
  fusions$transcript_id1 <- ifelse(rep(!("CDS_LEFT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_LEFT_ID)
  fusions$transcript_id2 <- ifelse(rep(!("CDS_RIGHT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_RIGHT_ID)
  fusions$fusion_transcript <- ifelse(rep(!("FUSION_CDS" %in% colnames(fusions)), nrow(fusions)), ".", toupper(sub("([a-z]*)", "\\1|", fusions$FUSION_CDS, perl=T)))
  fusions$reading_frame <- ifelse(rep(!("PROT_FUSION_TYPE" %in% colnames(fusions)), nrow(fusions)), ".", ifelse(fusions$PROT_FUSION_TYPE == "INFRAME", "in-frame", ifelse(fusions$PROT_FUSION_TYPE == "FRAMESHIFT", "out-of-frame", ".")))
  fusions$split_reads <- fusions$JunctionReadCount
  fusions$discordant_mates <- fusions$SpanningFragCount
  fusions$site1 <- rep("exon", nrow(fusions))
  fusions$site2 <- rep("exon", nrow(fusions))
  fusions$confidence <- rep("high", nrow(fusions))
  fusions$type <- ifelse(fusions$contig1 != fusions$contig2, "translocation", ifelse(fusions$direction1 == fusions$direction2, "inversion", ifelse((fusions$direction1 == "downstream") == (fusions$breakpoint1 < fusions$breakpoint2), "deletion", "duplication")))
} else {
  stop("Unrecognized fusion file format")
}

#pdf(outputFile, onefile=T, width=pdfWidth, height=pdfHeight, title=ifelse(sampleName != "", sampleName, fusionsFile))
#par(family=fontFamily)

if (nrow(fusions) == 0) {
  plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
  text(0, 0, "empty input file")
  warning("empty input file")
  dev.off()
  quit("no")
}


# insert dummy annotations for intergenic breakpoints:
#Agrega en el df de exons las zonas que son integenic:
if (any(fusions$site1 == "intergenic" | fusions$site2 == "intergenic")) {
  intergenicBreakpoints <- rbind(
    setNames(fusions[fusions$site1 == "intergenic",c("gene1", "strand1", "contig1", "breakpoint1")], c("gene", "strand", "contig", "breakpoint")),
    setNames(fusions[fusions$site2 == "intergenic",c("gene2", "strand2", "contig2", "breakpoint2")], c("gene", "strand", "contig", "breakpoint"))
  )
  exons <- rbind(exons, data.frame(
    contig=intergenicBreakpoints$contig,
    type="intergenic",
    start=sapply(intergenicBreakpoints$breakpoint-1000, max, 1),
    end=intergenicBreakpoints$breakpoint+1000,
    strand=".",
    attributes="",
    geneName=intergenicBreakpoints$gene,
    geneID=paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
    transcript=paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
    exonNumber="intergenic"
  ))
  fusions[fusions$site1 == "intergenic","gene_id1"] <- paste0(fusions[fusions$site1 == "intergenic","contig1"], ":", fusions[fusions$site1 == "intergenic","breakpoint1"])
  fusions[fusions$site2 == "intergenic","gene_id2"] <- paste0(fusions[fusions$site2 == "intergenic","contig2"], ":", fusions[fusions$site2 == "intergenic","breakpoint2"])
}



# layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
par(mar=c(0, 0, 0, 0))
#plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.4, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

#drawCoverage(left = min(exons1$left), right = max(exons1$right),
#             y = yCoverage, coverage = coverage1,
#             start = min(exons1$start), end = max(exons1$end), color = color1)

drawCoverage <- function(left, right, y, coverage, start, end, color) {

  maxResolution <- 5000 # max number of data points to draw coverage
  # draw coverage as bars
  if (!is.null(coverage)) {
    coverageData <- as.numeric(coverage[IRanges(sapply(start, max, min(start(coverage))), sapply(end, min, max(end(coverage))))])
    # downsample to maxResolution, if there are too many data points
    coverageData <- aggregate(coverageData, by=list(round(1:length(coverageData) * (right-left) * maxResolution/length(coverageData))), mean)$x
    #polygon(c(left, seq(left, right, length.out=length(coverageData)), right), c(y, y+coverageData*0.1, y), col=color, border=NA)
    x_vals <- seq(left, right, length.out = length(coverageData))
    y_vals <- y + coverageData * 0.3
    lines(x_vals, y_vals, col = color, lwd = 1)  # podés ajustar lwd (grosor)

  }
}

#drawStrand(left = breakpoint1, right = max(exonsOfGene$right), y = yFusion, color=darkColor1, strand = chartr("+-", "-+", exonsOfGene$strand[1]))
layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
par(mar=c(0, 0, 0, 0))
plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

drawStrand <- function(left, right, y, color, strand) {
  if (strand %in% c("+", "-")) {
    # draw strand
    lines(c(left+0.001, right-0.001), c(y, y), col=color, lwd=2)
    #lines(c(left+0.001, right-0.001), c(y, y), col=rgb(1,1,1,0.1), lwd=1)
    # indicate orientation
    if (right - left > 0.01)
      for (i in seq(left+0.005, right-0.005, by=sign(right-left-2*0.005)*0.05)) {
        arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=color, length=0.05, lwd=2, angle=60)
        arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=rgb(1,1,1,0.1), length=0.05, lwd=1, angle=60)
      }
  }
}

# layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
#layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
#par(mar=c(0, 0, 0, 0))
#plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.25, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

#drawExon(left = exons1[exon,"left"], right = exons1[exon,"right"], y = yExons, color = color1, title = exons1[exon,"exonNumber"], type = exons1[exon,"type"])

drawExon <- function(left, right, y, color, title, type) {
  exonHeight <- 0.03
  if (type == "CDS") {
    # draw coding regions as thicker bars
    rect(left, y+exonHeight, right, y+exonHeight/2-0.001, col=color, border=NA)
    rect(left, y-exonHeight, right, y-exonHeight/2+0.001, col=color, border=NA)
    # draw border
    lines(c(left, left, right, right), c(y+exonHeight/2, y+exonHeight, y+exonHeight, y+exonHeight/2), col=getDarkColor(color), lend=2)
    lines(c(left, left, right, right), c(y-exonHeight/2, y-exonHeight, y-exonHeight, y-exonHeight/2), col=getDarkColor(color), lend=2)

  } else if (type == "exon") {
    rect(left, y+exonHeight/2, right, y-exonHeight/2, col=color, border = NA)
    rect(left, y+exonHeight/2, right, y-exonHeight/2, col="white", border = getDarkColor(color))
    # add exon label
    text((left+right)/2, y, title, cex=0.9*fontSize, font = 1)
  }
}


drawProteinDomains <- function(fusion, exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors) {

  exonHeight <- 0.2
  exonsY <- 0.5
  geneNamesY <- exonsY - exonHeight/2 - 0.05

  # find coding exons
  codingExons1 <- exons1[exons1$type == "CDS" & fusion$site1 != "intergenic",]
  codingExons2 <- exons2[exons2$type == "CDS" & fusion$site2 != "intergenic",]

  # cut off coding regions beyond breakpoint
  if (fusion$direction1 == "upstream") {
    codingExons1 <- codingExons1[codingExons1$end >= fusion$breakpoint1,]
    codingExons1$start <- ifelse(codingExons1$start < fusion$breakpoint1, fusion$breakpoint1, codingExons1$start)
  } else {
    codingExons1 <- codingExons1[codingExons1$start <= fusion$breakpoint1,]
    codingExons1$end <- ifelse(codingExons1$end > fusion$breakpoint1, fusion$breakpoint1, codingExons1$end)
  }
  if (fusion$direction2 == "upstream") {
    codingExons2 <- codingExons2[codingExons2$end >= fusion$breakpoint2,]
    codingExons2$start <- ifelse(codingExons2$start < fusion$breakpoint2, fusion$breakpoint2, codingExons2$start)
  } else {
    codingExons2 <- codingExons2[codingExons2$start <= fusion$breakpoint2,]
    codingExons2$end <- ifelse(codingExons2$end > fusion$breakpoint2, fusion$breakpoint2, codingExons2$end)
  }

  # find overlapping domains
  exonsGRanges1 <- GRanges(codingExons1$contig, IRanges(codingExons1$start, codingExons1$end), strand=codingExons1$strand)
  exonsGRanges2 <- GRanges(codingExons2$contig, IRanges(codingExons2$start, codingExons2$end), strand=codingExons2$strand)
  domainsGRanges <- GRanges(proteinDomains$contig, IRanges(proteinDomains$start, proteinDomains$end), strand=proteinDomains$strand)
  domainsGRanges$proteinDomainName <- proteinDomains$proteinDomainName
  domainsGRanges$proteinDomainID <- proteinDomains$proteinDomainID
  domainsGRanges$color <- proteinDomains$color
  domainsGRanges <- domainsGRanges[suppressWarnings(unique(queryHits(findOverlaps(domainsGRanges, union(exonsGRanges1, exonsGRanges2)))))]

  # group overlapping domains by domain ID
  domainsGRangesList <- GRangesList(lapply(unique(domainsGRanges$proteinDomainID), function(x) { domainsGRanges[domainsGRanges$proteinDomainID == x] }))

  # trim protein domains to exon boundaries
  trimDomains <- function(domainsGRangesList, exonsGRanges) {
    do.call(
      "rbind",
      lapply(
        domainsGRangesList,
        function(x) {
          intersected <- as.data.frame(reduce(suppressWarnings(intersect(x, exonsGRanges))))
          if (nrow(intersected) > 0) {
            intersected$proteinDomainName <- head(x$proteinDomainName, 1)
            intersected$proteinDomainID <- head(x$proteinDomainID, 1)
            intersected$color <- head(x$color, 1)
          } else {
            intersected$proteinDomainName <- character()
            intersected$proteinDomainID <- character()
            intersected$color <- character()
          }
          return(intersected)
        }
      )
    )
  }
  retainedDomains1 <- trimDomains(domainsGRangesList, exonsGRanges1)
  retainedDomains2 <- trimDomains(domainsGRangesList, exonsGRanges2)

  # calculate length of coding exons
  codingExons1$length <- codingExons1$end - codingExons1$start + 1
  codingExons2$length <- codingExons2$end - codingExons2$start + 1

  # abort, if there are no coding regions
  if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
    text(0.5, 0.5, "Genes are not protein-coding.")
    return(NULL)
  }
  codingLength1 <- sum(codingExons1$length)
  codingLength2 <- sum(codingExons2$length)
  if (codingLength1 + codingLength2 == 0) {
    text(0.5, 0.5, "No coding regions retained in fusion transcript.")
    return(NULL)
  }
  if ((codingLength1 == 0 || grepl("\\.$", fusion$strand1)) && (codingLength2 == 0 || grepl("\\.$", fusion$strand2))) {
    text(0.5, 0.5, "Failed to determine retained protein domains due to lack of strand information.")
    return(NULL)
  }
  antisenseTranscription1 <- sub("/.*", "", fusion$strand1) != sub(".*/", "", fusion$strand1)
  antisenseTranscription2 <- sub("/.*", "", fusion$strand2) != sub(".*/", "", fusion$strand2)
  if ((codingLength1 == 0 || antisenseTranscription1) && (codingLength2 == 0 || antisenseTranscription2)) {
    text(0.5, 0.5, "No coding regions due to antisense transcription.")
    return(NULL)
  }

  # remove introns from protein domains
  removeIntronsFromProteinDomains <- function(codingExons, retainedDomains) {
    if (nrow(codingExons) == 0) return(NULL)
    cumulativeIntronLength <- 0
    previousExonEnd <- 0
    for (exon in 1:nrow(codingExons)) {
      if (codingExons[exon,"start"] > previousExonEnd)
        cumulativeIntronLength <- cumulativeIntronLength + codingExons[exon,"start"] - previousExonEnd
      domainsInExon <- which(between(retainedDomains$start, codingExons[exon,"start"], codingExons[exon,"end"]))
      retainedDomains[domainsInExon,"start"] <- retainedDomains[domainsInExon,"start"] - cumulativeIntronLength
      domainsInExon <- which(between(retainedDomains$end, codingExons[exon,"start"], codingExons[exon,"end"]))
      retainedDomains[domainsInExon,"end"] <- retainedDomains[domainsInExon,"end"] - cumulativeIntronLength
      previousExonEnd <- codingExons[exon,"end"]
    }
    # merge adjacent domains
    retainedDomains <- do.call(
      "rbind",
      lapply(
        unique(retainedDomains$proteinDomainID),
        function(x) {
          domain <- retainedDomains[retainedDomains$proteinDomainID == x,]
          merged <- reduce(GRanges(domain$seqnames, IRanges(domain$start, domain$end), strand=domain$strand))
          merged$proteinDomainName <- head(domain$proteinDomainName, 1)
          merged$proteinDomainID <- head(domain$proteinDomainID, 1)
          merged$color <- head(domain$color, 1)
          return(as.data.frame(merged))
        }
      )
    )
    return(retainedDomains)
  }
  retainedDomains1 <- removeIntronsFromProteinDomains(codingExons1, retainedDomains1)
  retainedDomains2 <- removeIntronsFromProteinDomains(codingExons2, retainedDomains2)

  # abort, if no domains are retained
  if (is.null(retainedDomains1) && is.null(retainedDomains2)) {
    text(0.5, 0.5, "No protein domains retained in fusion.")
    return(NULL)
  }

  # merge domains with similar coordinates
  mergeSimilarDomains <- function(domains, mergeDomainsOverlappingBy) {
    if (is.null(domains)) return(domains)
    merged <- domains[F,] # create empty data frame
    domains <- domains[order(domains$end - domains$start, decreasing=T),] # start with bigger domains => bigger domains are retained
    for (domain in rownames(domains)) {
      if (!any((abs(merged$start - domains[domain,"start"]) + abs(merged$end - domains[domain,"end"])) / (domains[domain,"end"] - domains[domain,"start"] + 1) <= 1-mergeDomainsOverlappingBy))
        merged <- rbind(merged, domains[domain,])
    }
    return(merged)
  }
  retainedDomains1 <- mergeSimilarDomains(retainedDomains1, mergeDomainsOverlappingBy)
  retainedDomains2 <- mergeSimilarDomains(retainedDomains2, mergeDomainsOverlappingBy)

  # if desired, reassign colors to protein domains to maximize contrast
  if (optimizeDomainColors) {
    uniqueDomains <- unique(c(retainedDomains1$proteinDomainID, retainedDomains2$proteinDomainID))
    # make rainbow of pretty pastell colors
    colors <- rainbow(length(uniqueDomains))
    colors <- apply(col2rgb(colors), 2, function(x) { 0.3 + y/255 * 0.7 }) # make pastell colors
    colors <- apply(colors, 2, function(x) {rgb(x["red"], x["green"], x["blue"])}) # convert back to rgb
    # reassign colors
    names(colors) <- uniqueDomains
    retainedDomains1$color <- colors[retainedDomains1$proteinDomainID]
    retainedDomains2$color <- colors[retainedDomains2$proteinDomainID]
  }

  # reverse exons and protein domains, if on the reverse strand
  if (any(codingExons1$strand == "-")) {
    codingExons1$length <- rev(codingExons1$length)
    temp <- retainedDomains1$end
    retainedDomains1$end <- codingLength1 - retainedDomains1$start
    retainedDomains1$start <- codingLength1 - temp
  }
  if (any(codingExons2$strand == "-")) {
    codingExons2$length <- rev(codingExons2$length)
    temp <- retainedDomains2$end
    retainedDomains2$end <- codingLength2 - retainedDomains2$start
    retainedDomains2$start <- codingLength2 - temp
  }

  # normalize length to 1
  codingExons1$length <- codingExons1$length / (codingLength1 + codingLength2)
  codingExons2$length <- codingExons2$length / (codingLength1 + codingLength2)
  retainedDomains1$start <- retainedDomains1$start / (codingLength1 + codingLength2)
  retainedDomains1$end <- retainedDomains1$end / (codingLength1 + codingLength2)
  retainedDomains2$start <- retainedDomains2$start / (codingLength1 + codingLength2)
  retainedDomains2$end <- retainedDomains2$end / (codingLength1 + codingLength2)

  # draw coding regions
  rect(0, exonsY-exonHeight/2, sum(codingExons1$length), exonsY+exonHeight/2, col=color1, border=NA)
  rect(sum(codingExons1$length), exonsY-exonHeight/2, sum(codingExons1$length) + sum(codingExons2$length), exonsY+exonHeight/2, col=color2, border=NA)

  # indicate exon boundaries as dotted lines
  exonBoundaries <- cumsum(c(codingExons1$length, codingExons2$length))
  if (length(exonBoundaries) > 1) {
    exonBoundaries <- exonBoundaries[1:(length(exonBoundaries)-1)]
    for (exonBoundary in exonBoundaries)
      lines(c(exonBoundary, exonBoundary), c(exonsY-exonHeight, exonsY+exonHeight), col="white", lty=3)
  }

  # find overlapping domains
  # nest if one is contained in another
  # stack if they overlap partially
  nestDomains <- function(domains) {
    if (length(unlist(domains)) == 0) return(domains)
    domains <- domains[order(domains$end - domains$start, decreasing=T),]
    rownames(domains) <- 1:nrow(domains)
    # find nested domains and make tree structure
    domains$parent <- 0
    for (domain in rownames(domains))
      domains[domains$start >= domains[domain,"start"] & domains$end <= domains[domain,"end"] & rownames(domains) != domain,"parent"] <- domain
    # find partially overlapping domains
    maxOverlappingDomains <- max(1, as.integer(coverage(IRanges(domains$start*10e6, domains$end*10e6))))
    padding <- 1 / maxOverlappingDomains * 0.4
    domains$y <- 0
    domains$height <- 0
    adjustPositionAndHeight <- function(parentDomain, y, height, padding, e) {
      for (domain in which(e$domains$parent == parentDomain)) {
        overlappingDomains <- which((between(e$domains$start, e$domains[domain,"start"], e$domains[domain,"end"]) |
                                       between(e$domains$end  , e$domains[domain,"start"], e$domains[domain,"end"])) &
                                      e$domains$parent == parentDomain)
        e$domains[domain,"height"] <- height/length(overlappingDomains) - padding * (length(overlappingDomains)-1) / length(overlappingDomains)
        e$domains[domain,"y"] <- y + (which(domain==overlappingDomains)-1) * (e$domains[domain,"height"] + padding)
        adjustPositionAndHeight(domain, e$domains[domain,"y"]+padding, e$domains[domain,"height"]-2*padding, padding, e)
      }
    }
    adjustPositionAndHeight(0, 0, 1, padding, environment())
    domains <- domains[order(domains$height, decreasing=T),] # draw nested domains last
    return(domains)
  }
  retainedDomains1 <- nestDomains(retainedDomains1)
  retainedDomains2 <- nestDomains(retainedDomains2)
  retainedDomains1$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains1$y
  retainedDomains2$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains2$y
  retainedDomains1$height <- retainedDomains1$height * (exonHeight-2*0.025)
  retainedDomains2$height <- retainedDomains2$height * (exonHeight-2*0.025)

  # draw domains
  drawProteinDomainRect <- function(left, bottom, right, top, color) {
    rect(left, bottom, right, top, col=color, border=getDarkColor(color))
    # draw gradients for 3D effect
    gradientSteps <- 20
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(top, bottom, len=gradientSteps), rgb(1,1,1,0.7))
    drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(bottom, bottom+(top-bottom)*0.4, len=gradientSteps), rgb(0,0,0,0.1))
  }
  if (length(unlist(retainedDomains1)) > 0)
    for (domain in 1:nrow(retainedDomains1))
      drawProteinDomainRect(retainedDomains1[domain,"start"], retainedDomains1[domain,"y"], retainedDomains1[domain,"end"], retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"], retainedDomains1[domain,"color"])
  if (length(unlist(retainedDomains2)) > 0)
    for (domain in 1:nrow(retainedDomains2))
      drawProteinDomainRect(sum(codingExons1$length)+retainedDomains2[domain,"start"], retainedDomains2[domain,"y"], sum(codingExons1$length)+retainedDomains2[domain,"end"], retainedDomains2[domain,"y"]+retainedDomains2[domain,"height"], retainedDomains2[domain,"color"])

  # draw gene names, if there are coding exons
  if (codingLength1 > 0)
    text(sum(codingExons1$length)/2, geneNamesY, fusion$gene1, font=2, cex=fontSize)
  if (codingLength2 > 0)
    text(sum(codingExons1$length)+sum(codingExons2$length)/2, geneNamesY, fusion$gene2, font=2, cex=fontSize)

  # calculate how many non-adjacent unique domains there are
  # we need this info to know where to place labels vertically
  countUniqueDomains <- function(domains) {
    uniqueDomains <- 0
    if (length(unlist(domains)) > 0) {
      uniqueDomains <- 1
      if (nrow(domains) > 1) {
        previousDomain <- domains[1,"proteinDomainID"]
        for (domain in 2:nrow(domains)) {
          if (previousDomain != domains[domain,"proteinDomainID"])
            uniqueDomains <- uniqueDomains + 1
          previousDomain <- domains[domain,"proteinDomainID"]
        }
      }
    }
    return(uniqueDomains)
  }
  if (length(unlist(retainedDomains1)) > 0)
    retainedDomains1 <- retainedDomains1[order(retainedDomains1$start),]
  uniqueDomains1 <- countUniqueDomains(retainedDomains1)
  if (length(unlist(retainedDomains2)) > 0)
    retainedDomains2 <- retainedDomains2[order(retainedDomains2$end, decreasing=T),]
  uniqueDomains2 <- countUniqueDomains(retainedDomains2)

  # draw title of plot
  titleY <- exonsY + exonHeight/2 + (uniqueDomains1 + 2) * 0.05
  text(0.5, titleY+0.01, "RETAINED PROTEIN DOMAINS", adj=c(0.5, 0), font=2, cex=fontSize)
  text(0.5, titleY, ifelse(fusion$reading_frame %in% c("in-frame", "out-of-frame"), paste(fusion$reading_frame, "fusion"), ifelse(fusion$reading_frame == "stop-codon", "stop codon before fusion junction", "reading frame unclear")), adj=c(0.5, 1), cex=fontSize)

  # draw domain labels for gene1
  if (length(unlist(retainedDomains1)) > 0) {
    previousConnectorX <- -1
    previousLabelX <- -1
    labelY <- exonsY + exonHeight/2 + uniqueDomains1 * 0.05
    for (domain in 1:nrow(retainedDomains1)) {
      # if possible avoid overlapping lines of labels
      connectorX <- min(retainedDomains1[domain,"start"] + 0.01, (retainedDomains1[domain,"start"] + retainedDomains1[domain,"end"])/2)
      if (connectorX - previousConnectorX < 0.01 && retainedDomains1[domain,"end"] > previousConnectorX + 0.01)
        connectorX <- previousConnectorX + 0.01
      labelX <- max(connectorX, previousLabelX) + 0.02
      # use a signle label for adjacent domains of same type
      adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains1) && retainedDomains1[domain+1,"proteinDomainID"] == retainedDomains1[domain,"proteinDomainID"]
      if (adjacentDomainsOfSameType) {
        labelX <- retainedDomains1[domain+1,"start"] + 0.015
      } else {
        text(labelX, labelY, retainedDomains1[domain,"proteinDomainName"], adj=c(0,0.5), col=getDarkColor(retainedDomains1[domain,"color"]), cex=fontSize)
      }
      lines(c(labelX-0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"]), col=getDarkColor(retainedDomains1[domain,"color"]))
      if (!adjacentDomainsOfSameType)
        labelY <- labelY - 0.05
      previousConnectorX <- connectorX
      previousLabelX <- labelX
    }
  }

  # draw domain labels for gene2
  if (length(unlist(retainedDomains2)) > 0) {
    previousConnectorX <- 100
    previousLabelX <- 100
    labelY <- exonsY - exonHeight/2 - (uniqueDomains2+1) * 0.05
    for (domain in 1:nrow(retainedDomains2)) {
      # if possible avoid overlapping connector lines of labels
      connectorX <- sum(codingExons1$length) + max(retainedDomains2[domain,"end"] - 0.01, (retainedDomains2[domain,"start"] + retainedDomains2[domain,"end"])/2)
      if (previousConnectorX - connectorX < 0.01 && sum(codingExons1$length) + retainedDomains2[domain,"start"] < previousConnectorX - 0.01)
        connectorX <- previousConnectorX - 0.01
      labelX <- min(connectorX, previousLabelX) - 0.02
      # use a signle label for adjacent domains of same type
      adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains2) && retainedDomains2[domain+1,"proteinDomainID"] == retainedDomains2[domain,"proteinDomainID"]
      if (adjacentDomainsOfSameType) {
        labelX <- sum(codingExons1$length) + retainedDomains2[domain+1,"end"] - 0.015
      } else {
        text(labelX, labelY, retainedDomains2[domain,"proteinDomainName"], adj=c(1,0.5), col=getDarkColor(retainedDomains2[domain,"color"]), cex=fontSize)
      }
      lines(c(labelX+0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains2[domain,"y"]), col=getDarkColor(retainedDomains2[domain,"color"]))
      if (!adjacentDomainsOfSameType)
        labelY <- labelY + 0.05
      previousConnectorX <- connectorX
      previousLabelX <- labelX
    }
  }

}

findExons <- function(exons, contig, geneID, direction, breakpoint, coverage, transcriptId, transcriptSelection) {
  # use the provided transcript if desired
  if (transcriptSelection == "provided" && transcriptId != "." && transcriptId != "") {
    candidateExons <- exons[exons$transcript == transcriptId,]
    if (nrow(candidateExons) == 0) {
      warning(paste0("Unknown transcript given in fusions file (", transcriptId, "), selecting a different one"))
    } else {
      return(candidateExons)
    }
  }

  if (transcriptSelection == "canonical") {
    candidateExons <- exons[exons$geneID == geneID & exons$contig == contig,]
  } else {
    # look for exon with breakpoint as splice site
    transcripts <- exons[exons$geneID == geneID & exons$contig == contig & exons$type == "exon" & (direction == "downstream" & abs(exons$end - breakpoint) <= 2 | direction == "upstream" & abs(exons$start - breakpoint) <= 2),"transcript"]
    candidateExons <- exons[exons$transcript %in% transcripts,]
    # if none was found, use all exons of the gene closest to the breakpoint
    if (nrow(candidateExons) == 0)
      candidateExons <- exons[exons$geneID == geneID & exons$contig == contig,]
    # if we have coverage information, use the transcript with the highest coverage if there are multiple hits
    if (!is.null(coverage)) {
      highestCoverage <- -1
      transcriptWithHighestCoverage <- NULL
      lengthOfTranscriptWithHighestCoverage <- 0
      for (transcript in unique(candidateExons$transcript)) {
        exonsOfTranscript <- candidateExons[candidateExons$transcript==transcript,]
        exonsOfTranscript$start <- sapply(exonsOfTranscript$start, max, min(start(coverage)))
        exonsOfTranscript$end <- sapply(exonsOfTranscript$end, min, max(end(coverage)))
        lengthOfTranscript <- sum(exonsOfTranscript$end - exonsOfTranscript$start + 1)
        coverageSum <- sum(as.numeric(coverage[IRanges(exonsOfTranscript$start, exonsOfTranscript$end)]))
        # we prefer shorter transcripts over longer ones, because otherwise there is a bias towards transcripts with long UTRs
        # => a longer transcript must have substantially higher coverage to replace a shorter one
        substantialDifference <- (1 - min(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage) / max(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage)) / 10
        if (lengthOfTranscript > lengthOfTranscriptWithHighestCoverage && coverageSum * (1-substantialDifference) > highestCoverage ||
            lengthOfTranscript < lengthOfTranscriptWithHighestCoverage && coverageSum > highestCoverage * (1-substantialDifference)) {
          highestCoverage <- coverageSum
          transcriptWithHighestCoverage <- transcript
          lengthOfTranscriptWithHighestCoverage <- lengthOfTranscript
        }
      }
      if (highestCoverage > 0)
        candidateExons <- candidateExons[candidateExons$transcript==transcriptWithHighestCoverage,]
    }
    # if the gene has multiple transcripts, search for transcripts which encompass the breakpoint
    if (length(unique(candidateExons$transcript)) > 1) {
      transcriptStart <- aggregate(candidateExons$start, by=list(candidateExons$transcript), min)
      rownames(transcriptStart) <- transcriptStart[,1]
      transcriptEnd <- aggregate(candidateExons$end, by=list(candidateExons$transcript), max)
      rownames(transcriptEnd) <- transcriptEnd[,1]
      encompassingExons <- between(breakpoint, transcriptStart[candidateExons$transcript,2], transcriptEnd[candidateExons$transcript,2])
      if (any(encompassingExons))
        candidateExons <- candidateExons[encompassingExons,]
    }
  }

  # find the consensus transcript, if there are multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    consensusTranscript <-
      ifelse(grepl("appris_principal_1", candidateExons$attributes), 12,
             ifelse(grepl("appris_principal_2", candidateExons$attributes), 11,
                    ifelse(grepl("appris_principal_3", candidateExons$attributes), 10,
                           ifelse(grepl("appris_principal_4", candidateExons$attributes), 9,
                                  ifelse(grepl("appris_principal_5", candidateExons$attributes), 8,
                                         ifelse(grepl("appris_principal", candidateExons$attributes), 7,
                                                ifelse(grepl("appris_candidate_longest", candidateExons$attributes), 6,
                                                       ifelse(grepl("appris_candidate", candidateExons$attributes), 5,
                                                              ifelse(grepl("appris_alternative_1", candidateExons$attributes), 4,
                                                                     ifelse(grepl("appris_alternative_2", candidateExons$attributes), 3,
                                                                            ifelse(grepl("appris_alternative", candidateExons$attributes), 2,
                                                                                   ifelse(grepl("CCDS", candidateExons$attributes), 1,
                                                                                          0
                                                                                   ))))))))))))
    candidateExons <- candidateExons[consensusTranscript == max(consensusTranscript),]
  }
  # use the transcript with the longest coding sequence, if there are still multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    codingSequenceLength <- ifelse(candidateExons$type == "CDS", candidateExons$end - candidateExons$start, 0)
    totalCodingSequenceLength <- aggregate(codingSequenceLength, by=list(candidateExons$transcript), sum)
    rownames(totalCodingSequenceLength) <- totalCodingSequenceLength[,1]
    candidateExons <- candidateExons[totalCodingSequenceLength[candidateExons$transcript,2] == max(totalCodingSequenceLength[,2]),]
  }
  # use the transcript with the longest overall sequence, if there are still multiple hits
  if (length(unique(candidateExons$transcript)) > 1) {
    exonLength <- candidateExons$end - candidateExons$start
    totalExonLength <- aggregate(exonLength, by=list(candidateExons$transcript), sum)
    rownames(totalExonLength) <- totalExonLength[,1]
    candidateExons <- candidateExons[totalExonLength[candidateExons$transcript,2] == max(totalExonLength[,2]),]
  }
  # if there are still multiple hits, select the first one
  candidateExons <- unique(candidateExons[candidateExons$transcript == head(unique(candidateExons$transcript), 1),])
  return(candidateExons)
}

findClosestGene <- function(exons, contig, breakpoint, extraConditions) {

  # find exons near breakpoint (extraConditions must define what is considered "near")
  closestExons <- exons[exons$contig == contig & extraConditions,] # find closest exon
  closestExons <- exons[exons$contig == contig & exons$geneID %in% closestExons$geneID,] # select all exons of closest gene

  # when more than one gene found with the given name, use the closest one
  if (length(unique(closestExons$geneID)) > 1) { # more than one gene found with the given name => use the closest one
    distanceToBreakpoint <- aggregate(1:nrow(closestExons), by=list(closestExons$geneID), function(x) { min(abs(closestExons[x,"start"]-breakpoint), abs(closestExons[x,"end"]-breakpoint)) })
    closestGene <- head(distanceToBreakpoint[distanceToBreakpoint[,2] == min(distanceToBreakpoint[,2]),1], 1)
    closestExons <- closestExons[closestExons$geneID == closestGene,]
  }

  # when no gene was found, return default values
  if (nrow(closestExons) == 0) {
    return(IRanges(max(1, breakpoint-1000), breakpoint+1000))
  } else {
    return(IRanges(min(closestExons$start), max(closestExons$end)))
  }
}

drawVerticalGradient <- function(left, right, y, color, selection=NULL) {
  # check if gradient should only be drawn in part of the region
  if (!is.null(selection)) {
    y <- y[selection]
    left <- left[selection]
    right <- right[selection]
  }
  # draw gradient
  for (i in 1:length(y)) {
    polygon(
      c(left[1:i], right[1:i]),
      c(y[1:i], y[i:1]),
      border=NA,
      col=rgb(col2rgb(color)["red",], col2rgb(color)["green",], col2rgb(color)["blue",], col2rgb(color, alpha=T)["alpha",]*(1/length(y)), max=255)
    )
  }
}
if (!render3dEffect) # nullify function, if no 3D effect should be drawn
  drawVerticalGradient <- function(left, right, y, color, selection=NULL) { }

drawCurlyBrace <- function(left, right, top, bottom, tip) {
  smoothness <- 20
  x <- cumsum(exp(-seq(-2.5, 2.5, len=smoothness)^2))
  x <- x/max(x)
  y <- seq(top, bottom, len=smoothness)
  lines(left+(tip-left)+x*(left-tip), y)
  lines(tip+x*(right-tip), y)
}

drawIdeogram <- function(adjust, left, right, y, cytobands, contig, breakpoint) {
  # define design of ideogram
  bandColors <- setNames(rgb(100:0, 100:0, 100:0, maxColorValue=100), paste0("gpos", 0:100))
  bandColors <- c(bandColors, gneg="#ffffff", acen="#ec4f4f", stalk="#0000ff")
  cytobands$color <- bandColors[cytobands$giemsa]
  arcSteps <- 30 # defines roundness of arc
  curlyBraceHeight <- 0.03
  ideogramHeight <- 0.04
  ideogramWidth <- 0.4
  # extract bands of given contig
  bands <- cytobands[cytobands$contig==contig,]
  if (nrow(bands) == 0) {
    warning(paste("Ideogram of contig", contig, "cannot be drawn, because no Giemsa staining information is available."))
    return(NULL)
  }
  # scale width of ideogram to fit inside given region
  bands$left <- bands$start / max(cytobands$end) * ideogramWidth
  bands$right <- bands$end / max(cytobands$end) * ideogramWidth
  # left/right-align cytobands
  offset <- ifelse(adjust=="left", left, right - max(bands$right))
  bands$left <- bands$left + offset
  bands$right <- bands$right + offset
  # draw curly braces
  tip <- min(bands$left) + (max(bands$right)-min(bands$left)) / (max(bands$end)-min(bands$start)) * breakpoint
  drawCurlyBrace(left, right, y-0.05+curlyBraceHeight, y-0.05, tip)
  # draw title of chromosome
  text((max(bands$right)+min(bands$left))/2, y+0.07, paste("chromosome", contig), font=2, cex=fontSize, adj=c(0.5,0))
  # draw name of band
  bandName <- bands[which(between(breakpoint, bands$start, bands$end)), "name"]
  text(tip, y+0.03, bandName, cex=fontSize, adj=c(0.5,0))
  # draw start of chromosome
  leftArcX <- bands[1,"left"] + (1+cos(seq(pi/2,1.5*pi,len=arcSteps))) * (bands[1,"right"]-bands[1,"left"])
  leftArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * (ideogramHeight/2)
  polygon(leftArcX, leftArcY, col=bands[1,"color"])
  # draw bands
  centromereStart <- NULL
  centromereEnd <- NULL
  for (band in 2:(nrow(bands)-1)) {
    if (bands[band,"giemsa"] != "acen") {
      rect(bands[band,"left"], y-ideogramHeight/2, bands[band,"right"], y+ideogramHeight/2, col=bands[band,"color"])
    } else { # draw centromere
      if (is.null(centromereStart)) {
        polygon(c(bands[band,"left"], bands[band,"right"], bands[band,"left"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
        centromereStart <- bands[band,"left"]
      } else {
        polygon(c(bands[band,"right"], bands[band,"left"], bands[band,"right"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
        centromereEnd <- bands[band,"right"]
      }
    }
  }
  # draw end of chromosome
  band <- nrow(bands)
  rightArcX <- bands[band,"right"] - (1+cos(seq(1.5*pi,pi/2,len=arcSteps))) * (bands[band,"right"]-bands[band,"left"])
  rightArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * ideogramHeight/2
  polygon(rightArcX, rightArcY, col=bands[band,"color"])
  # if there is no centromere, make an artificial one with length zero
  if (is.null(centromereStart) || is.null(centromereEnd)) {
    centromereStart <- bands[1,"right"]
    centromereEnd <- bands[1,"right"]
  }
  # draw gradients for 3D effect
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on p-arm
  drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on p-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on q-arm
  drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on q-arm
}


drawCircos <- function(fusion, fusions, cytobands, minConfidenceForCircosPlot, circosColors) {
  # check if Giemsa staining information is available
  for (contig in unlist(fusions[fusion,c("contig1", "contig2")])) {
    if (!any(cytobands$contig==contig)) {
      warning(paste0("Circos plot cannot be drawn, because no Giemsa staining information is available for contig ", contig, "."))
      # draw empty plots as placeholder
      plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
      plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
      return(NULL)
    }
  }
  # initialize with empty circos plot
  circos.clear()
  circos.initializeWithIdeogram(cytoband=cytobands, plotType=NULL)
  # use gene names as labels or <contig>:<position> for intergenic breakpoints
  geneLabels <- data.frame(
    contig=c(fusions[fusion,"contig1"], fusions[fusion,"contig2"]),
    start=c(fusions[fusion,"breakpoint1"], fusions[fusion,"breakpoint2"])
  )
  geneLabels$end <- geneLabels$start + 1
  geneLabels$gene <- c(fusions[fusion,"gene1"], fusions[fusion,"gene2"])
  geneLabels$gene <- ifelse(c(fusions[fusion,"site1"], fusions[fusion,"site2"]) == "intergenic", paste0(c(fusions[fusion,"display_contig1"], fusions[fusion,"display_contig2"]), ":", geneLabels$start), geneLabels$gene)
  # draw gene labels
  circos.genomicLabels(geneLabels, labels.column=4, side="outside", cex=fontSize, labels_height=0.27)
  # draw chromosome labels in connector plot
  for (contig in unique(cytobands$contig)) {
    set.current.cell(track.index=2, sector.index=contig) # draw in gene label connector track (track.index=2)
    circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex=0.85)
  }
  # draw ideograms
  circos.genomicIdeogram(cytoband=cytobands)
  # draw arcs
  confidenceRank <- c(low=0, medium=1, high=2)
  for (i in c(setdiff(1:nrow(fusions), fusion), fusion)) { # draw fusion of interest last, such that its arc is on top
    f <- fusions[i,]
    if (any(cytobands$contig==f$contig1) && any(cytobands$contig==f$contig2)) # ignore viral contigs, because we have no cytoband information for them
      if (minConfidenceForCircosPlot != "none" && confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot] || i==fusion)
        circos.link(
          f$contig1, f$breakpoint1,
          f$contig2, f$breakpoint2,
          lwd=2, col=ifelse(i==fusion, circosColors[f$type], getBrightColor(circosColors[f$type]))
        )
  }
  # draw legend
  plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
  legend(x="top", legend=names(circosColors), col=sapply(circosColors, getBrightColor), lwd=3, ncol=2, box.lty=0)
}

