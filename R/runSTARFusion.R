runSTAR <- function(patient_dir, twoPass = c("None","Basic")) {

  #STAR <- downloadSTAR()

  out <- downloadHG38()
  AnnotationHG38 <- out[[2]]
  FastaHG38 <- out[[1]]
  index_dir <- out[[3]]

  twoPass <- match.arg(twoPass[1], choices = c("None","Basic"))

  #if(missing(nThreads)){
  nThreads <- max(1,parallel::detectCores()-2)
  #}

  if(nThreads < 4){
    message(paste0("STAR running on ",nThreads, " may be not optimal"))
  }else{
    message(paste0("STAR running on ",nThreads))
  }


  file_list <- list.files(patient_dir)
  id <- basename(patient_dir)
  if(id  == "trimmed") {
    id <- basename(dirname(patient_dir))
  }

  #Evita repetir el analisis si ya fue hecho
  if (!(length(nchar(file_list[endsWith(file_list, sprintf("%s_Aligned_out.bam", id))])) == 0)) {
    if (is.na(file.info(file_list[endsWith(file_list, "Aligned_out.bam")])$size) ) {
      print("entre a eso de eliminar en la funcion del STAR")
      #file.remove(paste(patient_dir, file_list[endsWith(file_list, "Aligned_out.bam")], sep="/"))
      #file.remove(paste(patient_dir, file_list[endsWith(file_list, "Log.out")], sep ="/"))
      #file.remove(paste(patient_dir, file_list[endsWith(file_list, "Log.progress.out")] , sep= "/"))
    }
    message("The STAR alignment for this sample has already been done.")
    return(paste0(patient_dir, "/", file_list[endsWith(file_list, "Aligned_out.bam")], sep=""))
  }

  #Para que se pueda poner como entrada la carpeta de los trimmeados o la carpeta original:
  if  (basename(patient_dir) == "trimmed") {
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_1.fq.gz")], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "val_2.fq.gz")], sep="")
    patient_id <- basename(dirname(dirname(fileR1)))

  } else {
    fileR1 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "R1.fastq.gz")], sep="")
    fileR2 <- paste0(patient_dir, "/", file_list[endsWith(file_list, "R2.fastq.gz")], sep="")
    patient_id <- basename(dirname(fileR1))
  }

  if ((length(nchar(fileR1)) == 0) | (length(nchar(fileR2)) == 0)) {
    stop("There are no fastq files in this directory")
  }

  setwd(patient_dir)
  log <- system2(command = STAR,
                 args = c(paste0("--runThreadN ", nThreads),
                          paste0("--genomeDir " , index_dir),
                          paste0("--readFilesIn ", fileR1, " ", fileR2),
                          paste0("--readFilesCommand ", ifelse(stringr::str_detect(fileR1,".gz"),"zcat","-")),
                          # "--outStd BAM_Unsorted",
                          "--outSAMtype BAM Unsorted",
                          "--outSAMunmapped Within",
                          "--outBAMcompression 0",
                          "--outFilterMultimapNmax 50",
                          "--peOverlapNbasesMin 10",
                          "--alignSplicedMateMapLminOverLmate 0.5",
                          "--alignSJstitchMismatchNmax 5 -1 5 5",
                          "--chimSegmentMin 10",
                          "--chimOutType WithinBAM HardClip",
                          "--chimJunctionOverhangMin 10",
                          "--chimScoreDropMax 30",
                          "--chimScoreJunctionNonGTAG 0",
                          "--chimScoreSeparation 1",
                          "--chimSegmentReadGapMax 3",
                          "--chimMultimapNmax 50",
                          #sprintf("--outTmpDir %s/trimmed/%s_STARgenome", patient_dir, id ),
                          paste0("--outFileNamePrefix ", patient_id),
                          paste0("--twopassMode ", twoPass),
                          paste0("--sjdbGTFfile ", AnnotationHG38)
                 ), stdout = TRUE)#out.file)

  #Rename de bam file
  file.rename(sprintf("%s/%sAligned.out.bam", patient_dir, patient_id), sprintf("%s/%s_Aligned_out.bam", patient_dir, patient_id))

  message("STAR's analysis has finished!")

  return(sprintf("%s/%s_Aligned_out.bam", patient_dir, patient_id))
}





