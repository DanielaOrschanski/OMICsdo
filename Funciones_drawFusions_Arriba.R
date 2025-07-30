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
