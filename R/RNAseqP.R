#' @title RNAseqP
#' @description Executes the
#' @param patient_dir Path of the directory that contains R1 and R2 fastq files. It can be either the "trimmed" folder or the original folder. No tiene que tener espacios!!
#' @param genomeRef indicated which reference will be used. Write "HG38" or ...
#' @param fastQC_after_trim lfnrks
#' @param plot_FastQC_trim set to TRUE if you want to make a QC plot of the trimmed files.
#' @param R indicates if the curve in the QC plot will be constructed with R1 (R= "R1"), R2 (R= "R2") o both (R= "R1R2").
#' @param RunARRIBA set to TRUE if you want to run the ARRIBA software.
#' @param RunFeatureCounts vshbfj
#' @param RunMIXCR kvsjhbfj
#' @return softwares that were runned and the time they took.
#' @export
RNAseqP <- function(patients_dir,
                    genomeRef = "HG38",
                    fastQC_after_trim = FALSE,
                    plot_FastQC_trim = FALSE,
                    R = "R1R2",
                    trim_quality = 30,
                    RunARRIBA = FALSE,
                    RunFeatureCounts = TRUE,
                    RunMIXCR = FALSE
                    ) {

  #patients_dir <- "~/EnvironChile/Muestras"
  if (genomeRef == "HG38") {
    downloadHG38()
  }

  #Vectors that will be filled and show on a dataframe in the end of the function
  softwares_runned <- c()
  times_registered <- c()

  file_list <- list.dirs(path = patients_dir, full.names = FALSE, recursive = FALSE)

  #Functions applied to each patient one by one
  #patient <- file_list[1]
  for (patient in file_list) {
    patient_dir <- sprintf("%s/%s", patients_dir, patient)
    #FASTQC
    FastQC_time <- system.time(runFastQC(patient_dir))
    softwares_runned <- c(softwares_runned, "FastQC")
    times_registered <- c(times_registered, FastQC_time[[3]])
    message(sprintf("The patient %s has already been processed with FastQC", patient))

    #TRIMGALORE
    #patient_dir <- "~/EnvironChile/Muestras/Co1099"
      #Returns trimmed folder inside patient's folder
    TrimGalore_time <- system.time({
      patient_dir_trim <<- runTrimgalore(patient_dir, trim_quality = trim_quality)
    })
    times_registered <- c(times_registered, TrimGalore_time[[3]])
    softwares_runned <- c(softwares_runned, "TrimGalore")

    message(sprintf("The patient %s has already been processed with TrimGalore", patient))

    #FASTQC
    if (fastQC_after_trim == TRUE) {
      FastQC_trim_time <- system.time(runFastQC(patient_dir_trim))
      times_registered <- c(times_registered, FastQC_trim_time[[3]])
      softwares_runned <- c(softwares_runned, "FastQC after trim")
    }
    #patient_dir_trim <- "~/EnvironChile/Muestras/Co1109/trimmed"


    #MIXCR
    if(RunMIXCR == TRUE) {
      MIXCR_time <- system.time(runMIXCR(patient_dir_trim))
      times_registered <- c(times_registered, MIXCR_time[[3]])
      softwares_runned <- c(softwares_runned, "MIXCR")
      message(sprintf("The patient %s has already been processed with MIXCR", patient))
    }


    #STAR
    STAR_time <- system.time(runSTAR(patient_dir_trim))
    times_registered <- c(times_registered, STAR_time[[3]])
    softwares_runned <- c(softwares_runned, "STAR")
    message(sprintf("The patient %s has already been processed with STAR", patient))

    #ARRIBA
    if (RunARRIBA == TRUE) {
      tryCatch({
        if(genomeRef == "HG38") {
          genomeversion = "hg38"
          assemblyVersion = "GRCh38"
        } else {
          genomeversion = "hg19"
          assemblyVersion = "GRCh37"
        }
        ARRIBA_time <- runARRIBA(patient_dir_trim, genomeversion = genomeversion, assemblyVersion = assemblyVersion )
        #times_registered <- c(times_registered, ARRIBA_time[[3]])
        #softwares_runned <- c(softwares_runned, "ARRIBA")
        message(sprintf("The patient %s has already been processed with ARRIBA", patient))
      }, error = function(e) {
        message(sprintf("An error occurred while processing patient %s with ARRIBA: %s", patient, e$message))
      })
    }

  }

  #Functions applied to all the patients together:
  PLOT_time <- system.time(plotFastQC_PBSQ(patients_dir, trimmed = FALSE, R= R))
  times_registered <- c(times_registered, PLOT_time[[3]])
  softwares_runned <- c(softwares_runned, "plot_PBSQ")
  message("The PBSQ plot has been generated!")


  #Aca hay algo que da error:
  if(plot_FastQC_trim == TRUE) {
    PLOT_trim_time <- system.time(plotFastQC_PBSQ(patients_dir, trimmed = TRUE))
    times_registered <- c(times_registered, PLOT_trim_time[[3]])
    softwares_runned <- c(softwares_runned, "plot_PBSQ")
    message("The PBSQ plot has been generated with the trimmed fastqc!")
  }

  #FEATURE COUNTS
  if(RunFeatureCounts == TRUE) {
    #patients_dir <- "~/EnvironChile/Muestras"
    FC_time <- system.time({
      FC.object <- runFeatureCounts(patients_dir, genomeRef = genomeRef)
    })
    times_registered <- c(times_registered, FC_time)
    softwares_runned <- c(softwares_runned, "Feature_Counts")
    message("Feature Counts' analysis has finished!")
  }

  #TimeRegistration <- data.frame("Software" = softwares_runned, "Time" = times_registered)

  return(c(softwares_runned, times_registered))

}
