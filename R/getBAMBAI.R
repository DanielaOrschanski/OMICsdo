#' @title Get the BAM file from the patient selected by the path.
#' @description Check if there are BAM or BAI files already generated on the directory. If there are not, generates them.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @return bam.control path of the sorted BAM file
#' @export
getBAMBAI <- function(path_dir, ref = "HG19") {
  if (dir.exists(paste(path_dir, "/OMICsdo", sep = ""))) {
    file_list <- list.files(paste(path_dir, "/OMICsdo", sep=""))
    bam_file <- file_list[endsWith(file_list, "sortedR1R2.bam")]

    if (!(length(nchar(bam_file)) == 0)) {
      bam.control <- sprintf("%s/OMICsdo/%s", path_dir, bam_file)
    } else {
      bam.control <- fasta_to_bam(path_dir, ref = ref)
    }
    #bam.control <- "/media/16TBDisk/Daniela/CNVs/MuestrasCNVs/37336/OMICsdo/37336.recal.bam"
    bai_file <- file_list[endsWith(file_list, "sortedR1R2.bam.bai")]
    if (length(nchar(bai_file)) == 0) { indexBam(bam.control) }

  } else {
    bam.control <- fasta_to_bam(path_dir, ref = ref)
    indexBam(bam.control)
  }

  return(bam.control)
}


#' @title Generation of BAM file from FASTA or FASTQ files.
#' @description Given the path of a directory, this function will identify the R1 and R2 files
#' and will generate a SAM file, then a BAM and a sorted BAM file by applying commands from BWA,
#' Picard, Samtools and GATK and using the mitochondrial human genome reference.
#' All this softwares and reference are downloaded automathically with the installation of MitoR.
#' @param path_dir  path of the directory of the patient that will be analyzed.
#' @return path of the sorted BAM file.

fasta_to_bam <- function(path_dir, ref = "HG19") {

  omicsdo_sof <- sprintf("%s/OMICsdoSof", dirname(system.file(package = "OMICsdo")))
  softwares <- readLines(sprintf("%s/OMICsdoSof/path_to_soft.txt", dirname(system.file(package = "OMICsdo"))))

  #Checks all the softwares needed are downloaded ----
  linea_software <- grep("(?i)Samtools", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    Samtools <- downloadSamtools()
  } else {
    Samtools <-  sub("^Samtools ", "", linea_software)
  }

  linea_software <- grep("(?i)BWA", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    BWA <- downloadBWA()
  } else {
    BWA <-  sub("^BWA ", "", linea_software)
  }

  linea_software <- grep("(?i)PICARD", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    PICARD <- downloadPICARD()
  } else {
    PICARD <-  sub("^PICARD ", "", linea_software)
  }

  linea_software <- grep("(?i)GATK", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    GATK <- downloadGATK()
  } else {
    GATK <-  sub("^GATK ", "", linea_software)
  }

  if (ref == "HG38") {
    linea_software <- grep("(?i)HG38FASTA", softwares, ignore.case = TRUE, value = TRUE)
    if(length(nchar(linea_software)) == 0) {
      reference <- downloadHG38()
    } else {
      reference <- sub("^HG38FASTA ", "", linea_software)
    }
  } else {
    linea_software <- grep("(?i)HG19FASTA", softwares, ignore.case = TRUE, value = TRUE)
    if(length(nchar(linea_software)) == 0) {
      reference <- downloadHG19()
    } else {
      reference <- sub("^HG19FASTA ", "", linea_software)
    }
  }


  linea_software <- grep("(?i)TrimGalore", softwares, ignore.case = TRUE, value = TRUE)
  if(length(nchar(linea_software)) == 0) {
    TrimGalore <- downloadTrimGalore()
  } else {
    TrimGalore <- sub("^TrimGalore ", "", linea_software)
  }

  id <- basename(path_dir)
  file_list <- list.files(path_dir)
  patientR1 <- file_list[grep("R1.fastq.gz", file_list)]
  patientR2 <- file_list[grep("R2.fastq.gz", file_list)]
  patientR1 <- sprintf("%s/%s", path_dir, patientR1)
  patientR2 <- sprintf("%s/%s", path_dir,  patientR2)

  # Generate a folder inside the patient's folder where the files will be saved
  dir.create(sprintf("%s/OMICsdo", path_dir))

  #TrimGalore ----
  patient_dir_trim <- runTrimgalore(path_dir)

  #Alineamiento con BWA -----------------
  # Threads of the PC
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)
  cant_thread <- 15
  # Header-Contains information about the entire file
  #header <- sprintf('@RG\\tID:%s\\tLB:OMICsdo\\tSM:bar', id)

  # BWA mem - Alignment of R1 and R2 from the patient with the reference

  #Como paula
  #probar con mi referencia:
  #cd ~/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38
  #bwa index Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
  #cd /media/16TBDisk/Daniela/CNVs/MuestrasCNVs/37835
  #bwa mem /home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 37835_R1.fastq.gz 37835_R2.fastq.gz | gzip -3 > 37835_R1R2_miref.sam.gz


  #cd ~/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG19
  #bwa index Homo_sapiens.GRCh37.dna.primary_assembly.fasta

  patientR1 <- sprintf("%s/%s_val_1.fq.gz", patient_dir_trim, id)
  patientR2 <- sprintf("%s/%s_val_2.fq.gz", patient_dir_trim, id)

  #system2(BWA, sprintf("mem -t %s -o %s/OMICsdo/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, path_dir, id, reference, patientR1, patientR2, header), stdout = TRUE, wait = TRUE)

  #reference <- "/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
  system(sprintf("bwa mem -t %s %s %s %s  | gzip -3 > %s/OMICsdo/%s_R1R2.sam", cant_thread, reference, patientR1, patientR2, path_dir, id ))

  # SAM to BAM
  system2(Samtools, sprintf("view -bS %s/OMICsdo/%s_R1R2.sam > %s/OMICsdo/%s_R1R2.bam", path_dir, id, path_dir, id), stdout = TRUE, wait = TRUE)

  # SortSam (SortBAM)
  system2("java", sprintf("-jar %s SortSam -I %s/OMICsdo/%s_R1R2.bam -O %s/OMICsdo/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", PICARD, path_dir, id, path_dir, id), stdout = TRUE, wait = TRUE)

  #Remove big files
  #file.remove(sprintf("%s/OMICsdo/%s_R1R2.sam", path_dir, id))
  #file.remove(sprintf("%s/OMICsdo/%s_R1R2.bam", path_dir, id))

  return(sprintf("%s/OMICsdo/%s_sortedR1R2.bam", path_dir, id))
}
