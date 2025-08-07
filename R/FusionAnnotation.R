#ANOTACION DE DOMINIOS DE FUSIONES:
#library(OMICsdo)
#library(readr)
#Protein_domains <- read_tsv(file_Protein_domains_ARRIBA, comment = "#", col_names = FALSE)

#pfam_path <- downloadPfam(omicsdo_sof_path = "/media/16TBDisk/Daniela/OMICsdo/data")

#' @title downloadPfam
#' @description Generates
#' @param omicsdo_sof_path path
#' @export
#' @import openxlsx
#' @import readr
downloadPfam <- function(omicsdo_sof_path) {
  setwd(omicsdo_sof_path)

  pfam_path <- sprintf("%s/Pfam-A.hmm", omicsdo_sof_path)

  if(file.exists(pfam_path)) {
    message("Pfam was already downloaded")
  } else {
    setwd(dirname(pfam_path))
    message("Pfam database will be downloaded")
    system("wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz")
    system(sprintf("gunzip %s/Pfam-A.hmm.gz", dirname(pfam_path)))
    #system2("sudo apt install hmmer") # hacer esto en el terminal
    system(sprintf("hmmpress %s ", pfam_path))
  }
  return(pfam_path)
}

# Execution:
#protein_df <- readRDS("/media/respaldo8t/Daniela/Transcript_PeptideSequence.RDS")
#fusions_report <- read_excel("/media/8tb02/Daniela/Leucemia/Todos-FusionReports_Leucemia.xlsx")
#colnames(fusions_report)[1] <- "ID"

#fusions_report_anotated <- AnotateFusions(protein_df, fusions_report)
#fusions_report_anotated_domain <- GetDomains(fusions_report = fusions_report_anotated, path_dir = path_dir)
#domain = domain_1
#fusion = seq_prot


get_longest_common_substring <- function(domain, fusion) {
  domain_seq <- AAString(domain)
  fusion_seq <- AAString(fusion)

  hits <- pairwiseAlignment(domain_seq, fusion_seq, type = "local", scoreOnly = FALSE)
  conserved_domain <- as.character(subject(hits))
  return(conserved_domain)
}

#' @title GetDomains
#' @description Generates
#' @param fusions_report df
#' @param path_dir path to
#' @export
#' @import openxlsx
#' @import readr
GetDomains <- function(fusions_report, path_dir) {

  i=1
  Proteins_DB <- data.frame(Transcript_id = character(),Gene = character(), Sequence = character(), stringsAsFactors = FALSE)
  for (i in 1:nrow(fusions_report)) {

    # Añadir la primera fila con gene1 y sequence1
    Proteins_DB <- rbind(Proteins_DB, data.frame(Transcript_id = fusions_report$transcript_id1[i], Gene = fusions_report$gene1[i], Sequence = fusions_report$sequence1[i], stringsAsFactors = FALSE))
    # Añadir la segunda fila con gene2 y sequence2
    Proteins_DB <- rbind(Proteins_DB, data.frame(Transcript_id = fusions_report$transcript_id2[i], Gene = fusions_report$gene2[i], Sequence = fusions_report$sequence2[i], stringsAsFactors = FALSE))
  }

  Proteins_DB <- unique(Proteins_DB)
  if(any(Proteins_DB$Transcript_id == ".")) {
    Proteins_DB <- Proteins_DB[-which(Proteins_DB$Transcript_id == "."),]
  }
  if(any(Proteins_DB$Sequence == "-")) {
    Proteins_DB <- Proteins_DB[-which(Proteins_DB$Sequence == "-"),]
  }
  if(any(Proteins_DB$Sequence == "No se encuentra el transcript id1" | Proteins_DB$Sequence == "No se encuentra el transcript id2")) {
    Proteins_DB <- Proteins_DB[-which(Proteins_DB$Sequence == "No se encuentra el transcript id1" | Proteins_DB$Sequence == "No se encuentra el transcript id2"),]
  }

  #Escribo el fasta con todas esas proteinas:
  fasta_file <- file(sprintf("%s/fasta_seq_proteins", path_dir), open = "w")
  for (i in 1:nrow(Proteins_DB)) {
    # Escribir el Transcript_id y el Gene en la línea con el símbolo ">"
    writeLines(paste0(">", Proteins_DB$Transcript_id[i], "-", Proteins_DB$Gene[i]), fasta_file)

    # Escribir la secuencia correspondiente en la línea siguiente
    writeLines(Proteins_DB$Sequence[i], fasta_file)
  }

  close(fasta_file)

  # Scan domains with HMMR 3
  pfam_path <- downloadPfam(omicsdo_sof_path = "/media/16TBDisk/Daniela/OMICsdo/data")
  system(sprintf("hmmscan --domtblout %s/found-domains.tab %s %s/fasta_seq_proteins", path_dir, pfam_path, path_dir))
  system(sprintf("cat %s/found-domains.tab | grep -v '^#' | sed 's/  */\t/g' | cut -f 1,2,4,20,21 > %s/found-domains-extract.tab", path_dir, path_dir))

  domains_df <- read.table(sprintf("%s/found-domains-extract.tab", path_dir), header = FALSE, sep = "\t",
                           col.names = c("Domain", "DomainID", "Gen", "Start", "End"))

  domains_df <- domains_df[which(domains_df$Domain == "PK_Tyr_Ser-Thr"),]
  domains_df$TranscriptID <- sapply(strsplit(domains_df$Gen, split = "-"), `[`, 1)
  domains_df$Gen <- sapply(strsplit(domains_df$Gen, split = "-"), `[`, 2)

  i=1
  for(i in 1:nrow(domains_df)) {
    if(domains_df$TranscriptID[i] %in% fusions_report$transcript_id1) {
      domains_df$Sequence[i] <- fusions_report$sequence1[which(fusions_report$transcript_id1 == domains_df$TranscriptID[i])]
    } else if(domains_df$TranscriptID[i] %in% fusions_report$transcript_id2) {
      domains_df$Sequence[i] <- fusions_report$sequence2[which(fusions_report$transcript_id2 == domains_df$TranscriptID[i])]
    }

    domains_df$Domain_Sequence[i] <- substr(domains_df$Sequence[i], domains_df$Start[i], domains_df$End[i])

  }

  #Agregar coordenadas:
  #write.xlsx(domains_df, file = sprintf("%s/Domains_df.xlsx", path_dir))

  #Poner el dominio en el fusions_report: --------------------------------------
  i=1
  fusions_report$Domain1 <- "-"
  fusions_report$Domain2 <- "-"

  fusions_report$Perc_Conserved_Domain_1 <- "-"
  fusions_report$Conserved_Domain_1 <- "-"
  fusions_report$Coords_Conserved_Domain_1 <- "-"
  fusions_report$Residues_CBDOCK2_1 <- "-"

  fusions_report$Perc_Conserved_Domain_2 <- "-"
  fusions_report$Conserved_Domain_2 <- "-"
  fusions_report$Coords_Conserved_Domain_2 <- "-"
  fusions_report$Residues_CBDOCK2_2 <- "-"


  for( i in 1:nrow(fusions_report)) {
    seq_prot <- fusions_report$SeqProteinaFusion[i]

    #ANALIZO EL DOMINIO DEL GEN 1:
    if(fusions_report$transcript_id1[i] %in% domains_df$TranscriptID) {
      fusions_report$Domain1[i] <- domains_df$Domain_Sequence[which(domains_df$TranscriptID == fusions_report$transcript_id1[i])]
      domain_1 <- fusions_report$Domain1[i]

      if( grepl(domain_1, seq_prot) ) {
        message("El dominio 1 se encuentra al 100% en la secuencia de la proteina")
        porcentaje_conservacion_1 <- 100
        conserved_domain_1 <- domain_1

      } else {
        message("Analizando el porcentaje de conservación del dominio 1")
        library(Biostrings)
        conserved_domain_1 <- get_longest_common_substring(domain_1, seq_prot)
        #grepl(conserved_domain, fusions_report$SeqProteinaFusion[i])
        porcentaje_conservacion_1 <- round(nchar(conserved_domain_1) / nchar(domain_1) * 100)
      }

      fusions_report$Perc_Conserved_Domain_1[i] <- porcentaje_conservacion_1
      fusions_report$Conserved_Domain_1[i] <- conserved_domain_1

      # Encontrar coordenadas del conserved_dominio dentro de la proteina:
      start_domain <- regexpr(conserved_domain_1, seq_prot)[1]
      end_domain <- start_domain + nchar(conserved_domain_1) -1
      fusions_report$Coords_Conserved_Domain_1[i] <- paste0(start_domain, "-", end_domain)

      # Generar residuos para CBDOCK2:
      positions <- seq(start_domain, end_domain)
      amino_acids <- strsplit(domain_1, "")[[1]]
      fusions_report$Residues_CBDOCK2_1[i] <- paste0(positions, ":", amino_acids, collapse = ",")
    }


    #############################################################################
    if(fusions_report$transcript_id2[i] %in% domains_df$TranscriptID) { #----------------------------------------------

      fusions_report$Domain2[i] <- domains_df$Domain_Sequence[which(domains_df$TranscriptID == fusions_report$transcript_id2[i])]

      domain_2 <- fusions_report$Domain2[i]

      if( grepl(domain_2, seq_prot) ) {
        message("El dominio 2 se encuentra al 100% en la secuencia de la proteina")
        porcentaje_conservacion_2 <- 100
        conserved_domain_2 <- domain_2

      } else {
        message("Analizando el porcentaje de conservación del dominio 2")
        #Ver si por lo menos una parte esta conservadad  del dominio:
        library(Biostrings)
        conserved_domain_2 <- get_longest_common_substring(domain_2, seq_prot)
        #grepl(conserved_domain, fusions_report$SeqProteinaFusion[i])
        porcentaje_conservacion_2 <- round(nchar(conserved_domain_2) / nchar(domain_2) * 100)

      }

      fusions_report$Perc_Conserved_Domain_2[i] <- porcentaje_conservacion_2
      fusions_report$Conserved_Domain_2[i] <- conserved_domain_2

      # Encontrar coordenadas del conserved_dominio dentro de la proteina:
      start_domain <- regexpr(conserved_domain_2, seq_prot)[1]
      end_domain <- start_domain + nchar(conserved_domain_2) -1
      fusions_report$Coords_Conserved_Domain_2[i] <- paste0(start_domain, "-", end_domain)

      # Generar residuos para CBDOCK2:
      positions <- seq(start_domain, end_domain)
      amino_acids <- strsplit(domain_2, "")[[1]]
      fusions_report$Residues_CBDOCK2_2[i] <- paste0(positions, ":", amino_acids, collapse = ",")
    }

  }

  colnames(fusions_report)

  #fusions_report <- fusions_report[,c(1,2,3,4,5,6,7,8,9,10, 11,12,37, 38, 39, 13:36)]
  #write.xlsx(fusions_report, file = sprintf("%s/Fusions_Annot_Domains.xlsx", path_dir))
  return(fusions_report)
}

#ANOTACION DE FUSIONES -----------------------------------------------------------------------------------
#protein_df <- readRDS("/media/respaldo8t/Daniela/Transcript_PeptideSequence.RDS")
#fusions_report <- read_excel("/media/8tb02/Daniela/Leucemia/Todos-FusionReports_Leucemia.xlsx")
#colnames(fusions_report)[1] <- "ID"

#' @title AnotateFusions
#' @description Generates
#' @param fusions_report df
#' @param protein_df df
#' @export
#' @import readxl
#' @import readr
AnotateFusions <- function(protein_df, fusions_report) {

  fusions_report <- fusions_report[which(grepl("Protein_kinase_domain", fusions_report$retained_protein_domains)),]
  i=1
  fusions_report$SeqProteinaFusion <- "-"
  fusions_report$sequence1 <- "-"
  fusions_report$sequence2 <- "-"
  library(stringr)

  for (i in 1:nrow(fusions_report)) {

    transcript_id1 <- fusions_report$transcript_id1[i]
    agregar_AAfinal <- "NO"

    if(transcript_id1 %in% protein_df$Transcript) {
      fusions_report$sequence1[i] <- protein_df$Sequence[which(protein_df$Transcript == transcript_id1)]

      pedazo_seq1 <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][1]

      # funcion genera el predazo_seq_completo
      pedazo_seq1_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq1,
                                                         seq = fusions_report$sequence1[i],
                                                         gen = 1)

    } else { # SI NO ENCONTRAMOS EL EXACTO TRANSCRIPT ID: --------------------------

      fusions_report$sequence1[i] <- "No se encuentra el transcript id1"
      message("No se encuentra el transcript id1")
      transcript_id1_clean <- strsplit(transcript_id1, split="\\.")[[1]][1]

      if(transcript_id1_clean %in% protein_df$Transcript_clean) {
        message("SI se encuentra el transcript clean")
        nuevo_transcript_id1 <- protein_df$Transcript[which(protein_df$Transcript_clean == transcript_id1_clean)]

        fusions_report$sequence1[i] <- protein_df$Sequence[which(protein_df$Transcript == nuevo_transcript_id1)]
        pedazo_seq1 <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][1]

        # funcion genera el predazo_seq_completo
        pedazo_seq1_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq1,
                                                           seq = fusions_report$sequence1[i],
                                                           gen = 1)
      }
      else {
        message("NO se encuentra el transcript id1 ni en el clean")
        pedazo_seq1_completo <- "."
      }
    }

    #lo mismo para el trasncript id 2 ----------------------------------------
    transcript_id2 <- fusions_report$transcript_id2[i]

    if(transcript_id2 %in% protein_df$Transcript) {
      fusions_report$sequence2[i] <- protein_df$Sequence[which(protein_df$Transcript == transcript_id2)]
      agregar_AAfinal <- "NO"
      pedazo_seq2 <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][2]

      # funcion genera el predazo_seq_completo
      pedazo_seq2_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq2,
                                                         seq = fusions_report$sequence2[i],
                                                         gen = 2)

    } else {
      fusions_report$sequence2[i] <- "No se encuentra el transcript id2"
      message("No se encuentra el transcript id2")
      transcript_id2_clean <- strsplit(transcript_id2, split="\\.")[[1]][1]

      if(transcript_id2_clean %in% protein_df$Transcript_clean) {
        message("SI se encuentra el transcript clean")
        nuevo_transcript_id2 <- protein_df$Transcript[which(protein_df$Transcript_clean == transcript_id2_clean)]
        fusions_report$sequence2[i] <- protein_df$Sequence[which(protein_df$Transcript == nuevo_transcript_id2)]
        pedazo_seq2 <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][2]

        # funcion genera el predazo_seq_completo
        pedazo_seq2_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq2,
                                                           seq = fusions_report$sequence2[i],
                                                           gen = 2)
      } else {
        message("NO se encuentra el transcript id2 ni en el clean")
        pedazo_seq2_completo <- "."
      }
    }

    seq_prot_fusionada <- paste0(pedazo_seq1_completo, pedazo_seq2_completo)
    seq_prot_fusionada <- toupper(seq_prot_fusionada)
    fusions_report$SeqProteinaFusion[i] <- seq_prot_fusionada
  }
  return(fusions_report)
}


generate_PedazoSeqCompleto <- function(pedazo_seq, seq, gen) {
  if(is.na(pedazo_seq)) {
    message("pedazo es NA")
    pedazo_seq <- ""
  }
  if(!grepl(pedazo_seq, seq)) {
    message("Convierto minusuclas")
    pedazo_seq <- toupper(pedazo_seq)
  }
  if(grepl("\\?", pedazo_seq)) {
    pedazo_seq <- gsub("\\?", "", pedazo_seq)
  }

  if(grepl("\\*", pedazo_seq)) {
    message("Se encontró un codón de STOP")
    pedazo_seq <- gsub("\\*", "", pedazo_seq)
    pedazo_seq_completo <- ""
  } else {
    if(!grepl(pedazo_seq, seq)) {
      message("Elimino el ultimo AA")
      pedazo_seq <- substr(pedazo_seq, 1, nchar(pedazo_seq) - 1)
      agregar_AAfinal <- "SI"
    }
    if(!grepl(pedazo_seq, seq)) {
      message("No se encuentra el pedazo de peptide en la secuencia completa")
      pedazo_seq_completo <-  ""
    } else {
      message("Si se encontró el pedazo en la secuencia")
      #Si es gen 1 se queda con la parte de adelante de la secuencia.
      # Si es gen 2 se queda con la parte de atras de la secuencia
      pedazo_seq_completo <- ifelse(gen == 1,
                                    str_split(seq, pedazo_seq)[[1]][1],
                                    str_split(seq, pedazo_seq)[[1]][2]
      )
    }
  }

  if(is.na(pedazo_seq_completo) | pedazo_seq== "") {
    pedazo_seq_completo <- ""
  }

  #Si es gen 1 se pega secuencia + pedazo.
  # Si es gen 2 se pega pedazo + secuencia.
  pedazo_seq_completo <- ifelse(gen == 1,
                                paste0(pedazo_seq_completo, pedazo_seq),
                                paste0(pedazo_seq, pedazo_seq_completo))

  return(pedazo_seq_completo)
}
