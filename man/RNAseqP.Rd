% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAseqP.R
\name{RNAseqP}
\alias{RNAseqP}
\title{RNAseqP}
\usage{
RNAseqP(
  patients_dir,
  genomeRef = "HG38",
  fastQC_after_trim = FALSE,
  plot_FastQC_trim = FALSE,
  plotFastQC_PBSQ = FALSE,
  R = "R1R2",
  trim_quality = 30,
  RunARRIBA = TRUE,
  RunFeatureCounts = FALSE,
  RunMIXCR = FALSE
)
}
\arguments{
\item{genomeRef}{indicated which reference will be used. Write "HG38" or ...}

\item{fastQC_after_trim}{lfnrks}

\item{plot_FastQC_trim}{set to TRUE if you want to make a QC plot of the trimmed files.}

\item{R}{indicates if the curve in the QC plot will be constructed with R1 (R= "R1"), R2 (R= "R2") o both (R= "R1R2").}

\item{RunARRIBA}{set to TRUE if you want to run the ARRIBA software.}

\item{RunFeatureCounts}{vshbfj}

\item{RunMIXCR}{kvsjhbfj}

\item{patient_dir}{Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder. No tiene que tener espacios!!}
}
\value{
softwares that were runned and the time they took.
}
\description{
Executes the
}
