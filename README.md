# OMICsdo
A set of algorithms for DNA and RNA seq processing and its application in precision medicine.

# Installation
## Dependencies:
library(devtools)

install.packages("BiocManager")
BiocManager::install(c("Biostrings", "IRanges", "Rsamtools", "GenomicRanges", "GenomicAlignments"))

install_github(repo = "https://github.com/vplagnol/ExomeDepth.git")
library(ExomeDepth)

BiocManager::install("GEOquery")
install_github("wjawaid/enrichR")

install_github(repo = "https://github.com/DanielaOrschanski/OMICsdo.git")


## Authors

- Elmer Andrés Fernández (PhD) - Original Idea - [Profile](https://www.researchgate.net/profile/Elmer-Fernandez-2) - [Fundación para el Progreso ed la Medicina](https://fpmlab.org.ar/)- [CONICET](https://www.conicet.gov.ar)
- Daniela Orschanski - Developer and Maintainer - [Fundación para el Progreso ed la Medicina](https://fpmlab.org.ar/)- [CONICET](https://www.conicet.gov.ar)
