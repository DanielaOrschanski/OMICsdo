install.packages("BiocManager")
BiocManager::install("rtracklayer")

setwd("~/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba")
db_arriba <- "~/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3"
cyto <-"~/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv"
cyto <- read.delim(cyto)

db <- read.delim(db_arriba)
db2 <- db[which(grepl("Protein kinase domain", db$Name.Olfactory.receptor.color..0000FF.gene_id.ENSG00000186092.gene_name.OR4F5.protein_domain_id.PF13853)),]
colnames(db2)[9] <- "INFO"
library(tidyr)
db3 <- db2 %>%
  separate("INFO", into = c("col1", "col2", "col3", "col4", "col5"), sep = ";")
unique(db3$col5)

##########################################################
library(OMICsdo)

patient_dir <- "/media/4tb2/Daniela/Fusiones/Fusiones/MuestrasSRA/SRR15852401/trimmed"
#MAEQTYSc|qcrpcrtcslifmhsrspsqps

fasta_file <- "/media/4tb2/Daniela/Fusiones/Homo_sapiens.GRCh38.pep.all.fa.gz"
library(Biostrings)
protein_sequences <- readAAStringSet(fasta_file)
# Extraer los nombres y secuencias
protein_names <- names(protein_sequences)  # Nombres (contienen ID y descripción)
protein_seqs <- as.character(protein_sequences)  # Secuencias de proteínas

# Separar el ID y la descripción
split_names <- strsplit(protein_names, " ", fixed = TRUE)
protein_ids <- sapply(split_names, `[`, 1)  # Extrae el primer elemento como ID
protein_descriptions <- sapply(split_names, function(x) paste(x[-1], collapse = " "))  # Descripción

# Crear el dataframe
protein_df <- data.frame(
  ID = protein_ids,
  Description = protein_descriptions,
  Sequence = protein_seqs,
  stringsAsFactors = FALSE
)

protein_df$Description[1]
protein_df <- protein_df %>%
  separate(Description,
           into = c("Part1", "Part2", "Part3", "Transcript", "Rest"),
           sep = " ",
           extra = "merge", fill = "right") %>%
  # Extraer solo el transcript ID de la columna Transcript
  mutate(Transcript = sapply(strsplit(Transcript, " "), `[`, 1))

rownames(protein_df) <- NULL
protein_df$Transcript <- gsub("transcript:", "", protein_df$Transcript)

protein_df$Transcript_clean <- gsub("\\..*", "", protein_df$Transcript)
write_rds(protein_df, file = "/media/16TBDisk/Daniela/OMICsdo/data/Transcript_PeptideSequence.RDS")

# ACA EMPIEZA ·#########################################################################################################

#obtener residuos:
secuencia <- "NDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNS"
rango <- 97:157
residuos <- paste(rango, ":", unlist(strsplit(secuencia, "")), sep = "", collapse = ", ")


protein_df <- readRDS("/media/16TBDisk/Daniela/OMICsdo/data/Transcript_PeptideSequence.RDS")

library(readr)
fusions_report <- read_csv("/media/4tb2/Daniela/Fusiones/BRAFfusions.csv")
fusions_report <- read_excel("Fusiones/MuestrasSRA/SRR15852401/trimmed/SRR15852401_FusionReport.xlsx")
fusions_report <- read_excel("/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_FusionReport.xlsx")
#fusions_report <- read_excel("/media/16TBDisk/Daniela/Environ/TodoMama/Br0258/trimmed/Br0258_FusionReport.xlsx")
#fusions_report <- fusions_report[which(fusions_report$reading_frame =="in-frame"),]

#Me quedo con lo completo hasta el punto donde el ARRIBA dice que se corta:
library(readxl)

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
    if(!grepl(pedazo_seq1, fusions_report$sequence1[i])) {
      pedazo_seq1 <- toupper(pedazo_seq1)
    }
    if(!grepl(pedazo_seq1, fusions_report$sequence1[i])) {
      message("Elimino ultimo AA")
      agregar_AAfinal <- "SI"
      pedazo_seq1 <- substr(pedazo_seq1, 1, nchar(pedazo_seq1) - 1)
    }
    if(!grepl(pedazo_seq1, fusions_report$sequence1[i])) {
      #pedazo_seq2 <- substr(pedazo_seq2, 1, nchar(pedazo_seq2) - 1)
      message("No se encuentra el pedazo de peptide 1 en la secuencia completa")
    }

    pedazo_seq1_completo <- str_split(fusions_report$sequence1[i], pedazo_seq1)[[1]]
    if( agregar_AAfinal == "SI") {
      # Extraer el primer elemento y el primer carácter del segundo elemento
      pedazo_seq1_completo <- paste0(pedazo_seq1_completo[1], pedazo_seq1, substr(pedazo_seq1_completo[2], 1, 1))

      #Esto es por si me quedo con la mutacion: (osea con W)
      #mut <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][1]
      #ultima_letra <- toupper(substr(mut, nchar(mut), nchar(mut)))
      #pedazo_seq1_completo <- paste0(pedazo_seq1_completo[1], pedazo_seq1, ultima_letra)

    } else {
      pedazo_seq1_completo <- str_split(fusions_report$sequence1[i], pedazo_seq1)[[1]][1]
      pedazo_seq1_completo <- paste0(pedazo_seq1_completo, pedazo_seq1)
    }

    if(is.na(pedazo_seq1_completo)) {
      pedazo_seq1_completo <- ""
    }


  } else { # SI NO ENCONTRAMOS EL EXACTO TRANSCRIPT ID: --------------------------

    fusions_report$sequence1[i] <- "No se encuentra el transcript id1"
    message("No se encuentra el transcript id1")
    transcript_id1_clean <- strsplit(transcript_id1, split="\\.")[[1]][1]
    if(transcript_id1_clean %in% protein_df$Transcript_clean){
      message("SI se encuentra el transcript clean")
      nuevo_transcript_id1 <- protein_df$Transcript[which(protein_df$Transcript_clean == transcript_id1_clean)]
      #fusions_report$Warning[i] <- paste0("No se encuentra el transcript id1. Pero si se encuentra: ", nuevo_transcript_id1)

      fusions_report$sequence1[i] <- protein_df$Sequence[which(protein_df$Transcript == nuevo_transcript_id1)]
      pedazo_seq1 <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][1]

      pedazo_seq1_completo <- str_split(fusions_report$sequence1[i], pedazo_seq1)[[1]][1]
      pedazo_seq1_completo <- paste0(pedazo_seq1_completo, pedazo_seq1)
    }
  }

  #lo mismo para el trasncript id 2 ----------------------------------------
  transcript_id2 <- fusions_report$transcript_id2[i]

  if(transcript_id2 %in% protein_df$Transcript) {
    fusions_report$sequence2[i] <- protein_df$Sequence[which(protein_df$Transcript == transcript_id2)]
    agregar_AAfinal <- "NO"
    pedazo_seq2 <- str_split(fusions_report$peptide_sequence[i], "\\|")[[1]][2]

    if(is.na(pedazo_seq2)) {
      pedazo_seq2 <- ""
    }
    if(!grepl(pedazo_seq2, fusions_report$sequence2[i])) {
      message("Convierto minusuclas")
      pedazo_seq2 <- toupper(pedazo_seq2)
    }
    if(grepl("\\*", pedazo_seq2)) {
      message("Elimino *")
      pedazo_seq2 <- gsub("\\*", "", pedazo_seq2)
    }
    if(!grepl(pedazo_seq2, fusions_report$sequence2[i])) {
      message("Elimino el ultimo AA")
      pedazo_seq2 <- substr(pedazo_seq2, 1, nchar(pedazo_seq2) - 1)
      agregar_AAfinal <- "SI"
    }
    if(!grepl(pedazo_seq2, fusions_report$sequence2[i])) {
      #pedazo_seq2 <- substr(pedazo_seq2, 1, nchar(pedazo_seq2) - 1)
      message("No se encuentra el pedazo de peptide 2 en la secuencia completa")
    }

    #pedazo_seq2_completo <- str_split(fusions_report$sequence2[i], pedazo_seq2)[[1]]

    pedazo_seq2_completo <- str_split(fusions_report$sequence2[i], pedazo_seq2)[[1]][2]

    if(is.na(pedazo_seq2_completo) | pedazo_seq2== "") {
      pedazo_seq2_completo <- ""
    }

    pedazo_seq2_completo <- paste0(pedazo_seq2, pedazo_seq2_completo)

  } else {
    fusions_report$sequence2[i] <- "No se encuentra el transcript id2"
  }

  seq_prot_fusionada <- paste0(pedazo_seq1_completo, pedazo_seq2_completo)
  fusions_report$SeqProteinaFusion[i] <- seq_prot_fusionada
}

fusions_report <- fusions_report[,c(1,2, 21, 22, 23, 24, 5, 6, 29, 36, 37, 38, 7,8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35)]
#Generar una base de datos con Gen - transcript id - Secuencia DE INTERES
library(openxlsx)
write.xlsx(fusions_report, file = "/media/4tb2/Daniela/Fusiones/BRAF_fusions.xlsx")



######################################################################
# DOMINIOS

#setwd("/media/4tb2/Daniela/Fusiones")

#system("wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz")
#system("gunzip /media/4tb2/Daniela/Fusiones/Pfam-A.hmm.gz")
#sudo apt install hmmer
#system("hmmpress /media/4tb2/Daniela/Fusiones/Pfam-A.hmm")

#Generar un archivo con todas las proteinas de interes:
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
fasta_file <- file("/media/4tb2/Daniela/Fusiones/fasta_in", open = "w")
for (i in 1:nrow(Proteins_DB)) {
  # Escribir el Transcript_id y el Gene en la línea con el símbolo ">"
  writeLines(paste0(">", Proteins_DB$Transcript_id[i], "-", Proteins_DB$Gene[i]), fasta_file)

  # Escribir la secuencia correspondiente en la línea siguiente
  writeLines(Proteins_DB$Sequence[i], fasta_file)
}

close(fasta_file)


# Scan domains with HMMR 3
system("hmmscan --domtblout found-domains.tab Pfam-A.hmm fasta_in")

system("cat found-domains.tab | grep -v '^#' | sed 's/  */\t/g' | cut -f 1,2,4,20,21 > found-domains-extract.tab")


domains_df <- read.table("/media/4tb2/Daniela/Fusiones/found-domains-extract.tab", header = FALSE, sep = "\t",
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


write.xlsx(domains_df, file = "/media/4tb2/Daniela/Fusiones/Domains_df_BRAF.xlsx")

#Poner el dominio en el fusions_report:
i=1
fusions_report$Domain <- "-"
fusions_report$Coords_Domain <- "-"
fusions_report$Residues_CBDOCK2 <- "-"

for( i in 1:nrow(fusions_report)) {
  if(fusions_report$transcript_id2[i] %in% domains_df$TranscriptID) {
    fusions_report$Domain[i] <- domains_df$Domain_Sequence[which(domains_df$TranscriptID == fusions_report$transcript_id2[i])]
    grepl(fusions_report$Domain[i], fusions_report$SeqProteinaFusion[i])

    # Encontrar coordenadas del dominio dentro de la proteina:
    seq_proteina <- fusions_report$SeqProteinaFusion[i]
    domain <- fusions_report$Domain[i]
    start_domain <- regexpr(domain, seq_proteina)[1]
    end_domain <- start_domain + nchar(domain) -1
    fusions_report$Coords_Domain[i] <- paste0(start_domain, "-", end_domain)

    # Generar residuos para CBDOCK2:
    positions <- seq(start_domain, end_domain)
    amino_acids <- strsplit(domain, "")[[1]]
    fusions_report$Residues_CBDOCK2[i] <- paste0(positions, ":", amino_acids, collapse = ",")
  }
}

colnames(fusions_report)
fusions_report <- fusions_report[,c(1,2,3,4,5,6,7,8,9,10, 11,12,37, 38, 39, 13:36)]
write.xlsx(fusions_report, file = "/media/4tb2/Daniela/Fusiones/BRAF_fusions.xlsx")


#-------------------------------------------------------------------------

file_path <- "/media/4tb2/Daniela/Fusiones/Pfam-A.hmm"
hmm_data <- readLines(file_path)
head(hmm_data, 10)
hmm_data <- read.table(file_path, header = FALSE, sep = "", fill = TRUE, nrows = 100)

############################################################
#Leyendo un pdb:
install.packages('bio3d')
library(bio3d)

MLANApdb <- read.pdb("/media/4tb2/Daniela/Fusiones/PDBs/model_relaxed(1).pdb")
# Filtrar carbonos alfa
MLANA_carbonos_alfa <- MLANApdb$atom[MLANApdb$atom$elety == "CA", ]

# Filtrar los carbonos alfa
carbonos_alfa <- pdb$atom[pdb$atom$elety == "CA", ]
residuos <- MLANA_carbonos_alfa$resid
residuos_una_letra <- aa321(residuos)
cat(residuos_una_letra)

MLANApdb$seqres


KIApdb <- read.pdb("/media/4tb2/Daniela/Fusiones/PDBs/KIAA1549-BRAF-3.pdb")
# Filtrar carbonos alfa
KIA_carbonos_alfa <- KIApdb$atom[KIApdb$atom$elety == "CA", ]

BCRpdb <- read.pdb("/media/4tb2/Daniela/Fusiones/PDBs/BCR-ABL_e13a2.pdb")
pdb <- read.pdb("/media/4tb2/Daniela/Fusiones/PDBs/QKI-NTRK2-BP7243.pdb")
pdb <- read.pdb("/media/4tb2/Daniela/Fusiones/PDBs/QKI-NTRK2-BP1040.pdb")

BCRxyz <- BCRpdb$xyz
atoms <- pdb$atom
carbonos_alfa <- pdb$atom[pdb$atom$elety == "CA", ]
residuos <- carbonos_alfa$resid
residuos_una_letra <- aa321(residuos)
cat(residuos_una_letra)

####################################################
#

cif <- read.cif("/media/4tb2/Daniela/Fusiones/AlphaFold3-zip/fold_qki_ntrk2_bp7243_model_0.cif")

# Leer el archivo .cif como texto
ruta_cif <- "/media/4tb2/Daniela/Fusiones/AlphaFold3-zip/fold_qki_ntrk2_bp7243_model_0.cif"
lineas <- readLines(ruta_cif)

lineas_numeral <- grep("^#$", lineas)
lineas_atom_site <- grep("_atom_site", lineas)

inicio_seccion <- max(lineas_atom_site) +1
final_seccion <- max(lineas_numeral) - 1

seccion_atom_site <- lineas[inicio_seccion:final_seccion]


#Convertir en pdb:
library(dplyr)

#Agrego 2 lineas de terminacion:
lineas_vacias <- data.frame(
  group_PDB = rep(NA, 2),
  id = rep(NA, 2),
  type_symbol = rep(NA, 2),
  label_atom_id = rep(NA, 2),
  label_alt_id = rep(NA, 2),
  label_comp_id = rep(NA, 2),
  label_asym_id = rep(NA, 2),
  label_entity_id = rep(NA, 2),
  label_seq_id = rep(NA, 2),
  pdbx_PDB_ins_code = rep(NA, 2),
  Cartn_x = rep(NA, 2),
  Cartn_y = rep(NA, 2),
  Cartn_z = rep(NA, 2),
  occupancy = rep(NA, 2),
  B_iso_or_equiv = rep(NA, 2),
  auth_seq_id = rep(NA, 2),
  auth_asym_id = rep(NA, 2),
  pdbx_PDB_model_num = rep(NA, 2)
)
df_atom_site_completo <- rbind(df_atom_site, lineas_vacias)
df_atom_site_completo$group_PDB[nrow(df_atom_site_completo)] <- "END"
df_atom_site_completo$group_PDB[nrow(df_atom_site_completo)-1] <- "TER"
df_atom_site_completo$id[nrow(df_atom_site_completo)-1] <- df_atom_site_completo$id[nrow(df_atom_site_completo)-2] +1
df_atom_site_completo$label_comp_id[nrow(df_atom_site_completo)-1] <-df_atom_site_completo$label_comp_id[nrow(df_atom_site_completo)-2]
df_atom_site_completo$label_asym_id[nrow(df_atom_site_completo)-1] <-df_atom_site_completo$label_asym_id[nrow(df_atom_site_completo)-2]


CA <- df_atom_site_completo[which(df_atom_site_completo$label_atom_id == "CA"),]
library(bio3d)

# Suponiendo que df_atom_site_completo ya está cargado y contiene las columnas mencionadas:
# Filtrar filas donde alguna de las coordenadas Cartn_x, Cartn_y o Cartn_z sea NA
df_atom_site_completo <- df_atom_site_completo[!is.na(df_atom_site_completo$Cartn_x) &
                                                 !is.na(df_atom_site_completo$Cartn_y) &
                                                 !is.na(df_atom_site_completo$Cartn_z), ]

xyz <- as.matrix(df_atom_site_completo[, c("Cartn_x", "Cartn_y", "Cartn_z")])  # Coordenadas xyz
colnames(xyz) <- c("x", "y", "z")
any(is.na(xyz))

type <- rep("ATOM", nrow(df_atom_site_completo))  # Tipo de registro (suponiendo que son todos átomos)
resno <- df_atom_site_completo$label_seq_id  # Número de residuo
resid <- df_atom_site_completo$label_comp_id  # Tipo de residuo
eleno <- df_atom_site_completo$id  # Número de átomo
elety <- df_atom_site_completo$type_symbol  # Tipo de átomo (por ejemplo, C, N, etc.)
chain <- df_atom_site_completo$label_asym_id  # Identificador de cadena
insert <- rep("", nrow(df_atom_site_completo))  # Insertion code (vacío si no se especifica)
alt <- rep("", nrow(df_atom_site_completo))  # Alternativa (vacío si no se especifica)
o <- rep(1.0, nrow(df_atom_site_completo))  # Ocupación
b <- df_atom_site_completo$B_iso_or_equiv  # B-factor
segid <- rep("", nrow(df_atom_site_completo))  # Segment ID (vacío si no se especifica)
elesy <- df_atom_site_completo$type_symbol  # Símbolo del elemento (lo mismo que 'elety')
charge <- rep(NA, nrow(df_atom_site_completo))  # Carga atómica (vacío si no se especifica)

# Crear el archivo PDB
write.pdb(
  pdb = NULL,
  file = "/media/4tb2/Daniela/Fusiones/PDBs/QKI_NTRK2_bp7243_desdeCIF.pdb",
  xyz = xyz,
  type = type,
  resno = resno,
  resid = resid,
  eleno = eleno,
  elety = elety,
  chain = chain,
  insert = insert,
  alt = alt,
  o = o,
  b = b,
  segid = segid,
  elesy = elesy,
  charge = charge,
  append = FALSE,
  verbose = TRUE,
  chainter = TRUE,
  end = TRUE
)



#writeLines(as.matrix(df_atom_site_completo), con = "/media/4tb2/Daniela/Fusiones/PDBs/QKI_NTRK2_bp7243_desdeCIF.pdb")
# Escribir el archivo PDB de manera más simple, con un formato tabular ajustado
write.table(df_atom_site_completo,
            file = "/media/4tb2/Daniela/Fusiones/PDBs/QKI_NTRK2_bp7243_desdeCIF.pdb",
            quote = FALSE,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE)



pdb <- read.pdb("/media/4tb2/Daniela/Fusiones/PDBs/fold_qki_ntrk2_bp7243_model_0.pdb")
pdb$atom
carbonos_alfa <- pdb$atom[pdb$atom$elety == "CA", ]
residuos <- MLANA_carbonos_alfa$resid
residuos_una_letra <- aa321(residuos)
cat(residuos_una_letra)
