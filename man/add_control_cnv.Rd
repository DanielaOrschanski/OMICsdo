% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_control_cnv.R
\name{add_control_cnv}
\alias{add_control_cnv}
\title{Add patient to Control DB}
\usage{
add_control_cnv(path_dir, bed, ref)
}
\arguments{
\item{path_dir}{path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta/q files).}
}
\value{
Control DB updated
}
\description{
Saves reads' counts of each gene from the sequence of the patient in a reference data set use as a control data base.
}
\examples{
ControlDB <- AddControl(path_dir ="/home/sam/Patients/123/123.fastq")
}
