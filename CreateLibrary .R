install.packages("miniUI", dependencies = TRUE)
install.packages("pkgdown", dependencies = TRUE)
#install.packages("systemfonts", type = "source", repos = "http://cran.us.r-project.org", INSTALL_opts = "--no-lock", dependencies = TRUE)

#sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-devfrd

install.packages("devtools", dependencies = TRUE)
library(devtools)
#devtools::create(path = "/home/daniela/OMICsdo")

path.expand(Sys.getenv("R_LIBS_USER"))
"/home/juan/R/x86_64-pc-linux-gnu-library/4.1"
setwd("/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof")

setwd("/media/16TBDisk/Daniela/OMICsdo")
setwd("/media/4tb2/Daniela/Biota/PipelineBiota")
devtools::build()
#devtools::load_all() #para asegurarte de que se está ejecutando correctamente
devtools::document()
setwd("/media/16TBDisk/Daniela/OMICsdo")
setwd("/media/4tb2/Daniela/Biota/PipelineBiota")
devtools::install()
#devtools::install_local("/home/daniela/OMICsdo_0.0.0.9000.tar.gz")

library(PipelineBiota)

#Errores en instalacion:
txt <- read_file("~/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/path_to_soft.txt")

unlink("/home/juan/R/x86_64-pc-linux-gnu-library/4.1/00LOCK-OMICsdo", recursive = TRUE)
unlink("/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdo", recursive = TRUE)
#restart R

library(OMICsdo)

patients_dir <- "/media/16TBDisk/Daniela/TodoMama"
patients_dir <- "/media/16TBDisk/Daniela/58MuestrasFaltantes/Mama"
RNAseqP(patients_dir,
        genomeRef = "HG38",
        fastQC_after_trim = TRUE,
        plot_FastQC_trim = TRUE,
        trim_quality = 20,
        RunARRIBA = FALSE,
        RunFeatureCounts = TRUE,
        RunMIXCR = FALSE
        )

runFastQC("/media/16TBDisk/Daniela/TodoMama/Br0587")

patient_dir <- "/media/16TBDisk/Daniela/Biota/MuestrasLeo/728"
runFastQC(patient_dir)

#Error con trimgalore: falta el cutadapt
#sudo apt install cutadapt

#Error en los nombres de los archivos de salida del STAR
#ERROR EN EL STAR con el index de la referencia
downloadHG38()
#Creo que era porque el path que se guardo de la annotacion y del fastq tenian el .gz al final

#STAR download and installation completed successfully
ESTOY EN DOWNLOAD HG38
The annotation for HG38 was already downloaded
The fasta file for HG38 was already downloaded
ESTOY EN EL INDEX
The reference genome was already indexed
[2024-01-16T11:58:12] Launching Arriba 2.4.0
ERROR: file not found/readable: /home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz

FUSION FILES NOT FOUND
/media/16TBDisk/Daniela/MuestrasColon/Co1100/trimmed/Co1100_fusions.tsv
/media/16TBDisk/Daniela/MuestrasColon/Co1100/trimmed/Co1100_fusions_discarded.tsvError in file(file, "rt") : no se puede abrir la conexión
Además: Warning message:
  In file(file, "rt") :
  no fue posible abrir el archivo '/media/16TBDisk/Daniela/MuestrasColon/Co1100/trimmed/Co1100_fusions.tsv': No such file or directory
Timing stopped at: 0.006 0 0.007





counts_viejas <- as.data.frame(FeatureCount_Report_56_mama$counts)
count_0258_v <- as.data.frame(counts_viejas$Br0258_Aligned_out.bam)


counts_nuevas <- as.data.frame(FeatureCount_Report$counts)
count_0258_n <- subset(counts_nuevas, select ="Br0258_Aligned_out.bam" )
save(count_0258_n, file = "/media/16TBDisk/Daniela/TodoMama/Counts_0258_nuevo.RData")
