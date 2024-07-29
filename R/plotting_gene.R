#' @title Generate a PileUp plot of the mutated genes and the mutations themselves.
#' @description Create a visual representation of the read coverage and amount of mutations per gene.
#' Because of the large amount of positions on each gene, it also creates a pileup plot for each mutation and a range of near nucleotides.
#' @param BAM_file Patient's BAM file.
#' @param plot_mutation Whether to plot each mutation or not
#' @param range Amount of nucleotides to show next to the mutated one. plot_mutation must be TRUE in this case.
#' @return Nested list with the plots. The main list is named after the mutated genes. Inside of them there is the gene plot and all mutation plots belonging to that gene.
#' @examples
#' pileupPlot("47286", plot_mutation = FALSE)
#' pileupPlot("47286", plot_mutation = TRUE)
#' pileupPlot("47286", plot_mutation = TRUE, range = 10)
#' @import showtext
#' @import ggtext

pileupPlot <- function(BAM_file, plot_mutation = TRUE, range = 5) {
  # Takes the Patient's ID from the BAM file path
  patient_ID <- basename(BAM_file) %>% strsplit("_") %>% unlist()
  patient_ID <- patient_ID[1]
  
  # Checks which genes have mutations. Result is a DF with colomns: Interval, Gene, Mutations
  mutations_at_genes <- generateTXT_gene(patient_ID, range)
  
  mitor_sof <- sprintf("%s/mitorDB/Softwares", Sys.getenv('R_LIBS_USER'))
  GATK <- sprintf("%s/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
  reference <- sprintf("%s/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_sof)
  
  path_to_MitoR <- substr(BAM_file, start = 0, stop = (nchar(BAM_file)-nchar(basename(BAM_file))))
  
  plots_to_xlsx <- list()
  showtext_auto()
  for (i in (1:nrow(mutations_at_genes))) {
    
    # GATK function to get the Allelic count of the listed positions. In this case, the complete mutated gene
    system2("java", sprintf("-jar %s CollectAllelicCounts -R %s -O %s%s_region_allelic -I %s -L %s", GATK, reference, path_to_MitoR, patient_ID, BAM_file, mutations_at_genes$Interval[i]), stdout = TRUE)
    
    # Takes the information given by GATK and creates a DF with it. Colomns: POS, REF Counts, ALT Counts, REF, ALT
    allelic_coverage <- readLines(sprintf("%s%s_region_allelic", path_to_MitoR, patient_ID))[c(-1,-2,-3, -4)] %>%
      strsplit("\t") %>% unlist()
    allelic_coverage <- matrix(allelic_coverage, nrow = 6, ncol = (length(allelic_coverage)/6)) %>% t() %>%
      as.data.frame()
    colnames(allelic_coverage) <- c("CHROM", "POS", "REF Counts", "ALT Counts", "REF", "ALT")
    allelic_coverage <- allelic_coverage[, c("POS", "REF Counts", "ALT Counts", "REF", "ALT")]
    
    # Cells with numbers to integers
    allelic_coverage$`REF Counts` <- as.integer(allelic_coverage$`REF Counts`)
    allelic_coverage$`ALT Counts` <- as.integer(allelic_coverage$`ALT Counts`)
    
    # Creates the colomn of Total counts (ALT Counts + REF Counts)
    allelic_coverage <- cbind(allelic_coverage, "Total Counts" = (allelic_coverage$`REF Counts` + allelic_coverage$`ALT Counts`))
    
    # Designs the gene plot and saves it on the list
    gene_plot <- plot_gene_design(allelic_coverage, mutations_at_genes[i, ])
    
    # The name of the element will be the gen name
    plots_to_xlsx[[mutations_at_genes$gene[i]]][[mutations_at_genes$gene[i]]] <- gene_plot
    
    # PLOTMUT---------
    if (plot_mutation == TRUE) {
      
      # Takes the mutations to plot from the DF
      mutations_to_plot <- mutations_at_genes$mut[i] %>%
        strsplit(" ") %>% unlist()
      
      for (mutation in mutations_to_plot) {
        
        # Takes the POS and ALT from the mutation to plot
        params <- strsplit(mutation, ">") %>% unlist()
        position <- params[2]
        alternative_allele <- params[3]
        
        # Extract the meaningful (depending on the range) rows of the DF for the mutation plot
        inferior_row <- which(allelic_coverage$POS == (as.integer(position) - range))
        superior_row <- which(allelic_coverage$POS == (as.integer(position) + range))

        data_mutation_plot <- allelic_coverage[inferior_row : superior_row, ]
        
        mutation_plot <- plot_mutation_design(data_mutation_plot, alternative_allele, range)
        
        # Each mutation will be saved inside the gen list element
        plots_to_xlsx[[mutations_at_genes$gene[i]]][[position]] <- mutation_plot
      }
    }
  }
  file.remove(sprintf("%s%s_dedup_reads.bam", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_dedup_reads.bai", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_region_allelic", path_to_MitoR, patient_ID))
  
  return(plots_to_xlsx)
}

#' @title Check which genes have mutations.
#' @description Using a universal BED file and the mutations of the patient detected by MitoR, the function gives a list
#' of mutations, which gene has them and the interval of those genes.
#' @param patient_ID Patient ID.
#' @return Updated BED table as a DataFrame object. Containing colomns: Interval, gene name, mutations found (REF-POS-ALT)
#' @import showtext

generateTXT_gene <- function(patient_ID, range) {
  # Load mutations
  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  path <- as.character(readRDS(sprintf("%s/RDS_DB.rds", mitor_db))[[sprintf("%s", patient_ID)]][1])
  mutations <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))[[sprintf("%s", patient_ID)]]$Mutations
  mutations <-  stringr::str_split(mutations,"/", simplify= TRUE) %>% as.data.frame(row.names = NULL)
  colnames(mutations) <- c("ALT", "POS", "REF")

  # Load BED file with: Gene name, start position, end position
  gene_start_end <- bedSNP
  
  # Dealing with insertions and deletions on the length of the gene
  gene_start_end <- indels_in_gene(mutations, gene_start_end)
  
  # Looks for the genes with mutations in them
  mutations_at_genes <- genes_of_mutations(mutations, gene_start_end, range)
  
  return(mutations_at_genes)
}


#' @title Modifies the BED table. Makes a gene longer or shorter depending on the patients variants.
#' @description Depending on the position of the insertion or deletion of the variant, it changes the length of the gene's start and end
#' @param mutations DataFrame of colomns REF-POS-ALT
#' @param start_end_gene BED file 
#' @return Updated BED table as a DataFrame object
indels_in_gene <- function(mutations, start_end_gene) {
  
  j <- 1
  for (i in (1:nrow(mutations))) {
    indel <- nchar(mutations[i, "ALT"]) - nchar(mutations[i, "REF"])
    if (indel == 0) {
      # Continues
    } else {
      # Insertion or Deletion      
      while(!((as.integer(mutations[i, "POS"]) >= as.integer(start_end_gene[j, "start"])) && (as.integer(mutations[i, "POS"]) <= as.integer(start_end_gene[j, "end"])))){
        j <- j + 1
      }
      start_end_gene[j, "end"] <- start_end_gene[j, "end"] + indel
      if (j < nrow(start_end_gene)){
        start_end_gene[j+1, "start"] <- start_end_gene[j+1, "start"] + indel
      }
    }
  }
  return(start_end_gene)
}

#' @title Modifies the patients mutation table.
#' @description Makes modifications to the mutation's table
#' @param mutations DataFrame of colomns REF-POS-ALT
#' @param gene_start_end BED file with modified length depending on the insertions/deletions of the patient
#' @return Updated BED table as a DataFrame object. Containing colomns: Interval, gene name, mutations found (REF-POS-ALT)

genes_of_mutations <- function(mutations, gene_start_end, range) {
  position_of_mutations <- unlist(as.integer(mutations$POS))

  i <- 1 # POS row
  j <- 1 # Gene Row
  genes_of_mutations <- data.frame()
  
  while (j < nrow(gene_start_end)) {
    while ((position_of_mutations[i] > gene_start_end$start[j]) && (position_of_mutations[i] < gene_start_end$end[j])) {
      if (gene_start_end$name[j] == "non-coding") {
        NULL
      } else if (!(gene_start_end$name[j] %in% genes_of_mutations$name)) {
        genes_of_mutations <- rbind(genes_of_mutations, cbind(gene_start_end[j, ], mut = paste0(mutations[i,], sep = "", collapse = ">")))
        # Add it. Don't go to the next row because it might not be the only mutation in this gene
      } else {
        genes_of_mutations[nrow(genes_of_mutations), "mut"] <- paste(genes_of_mutations[nrow(genes_of_mutations), "mut"], paste0(mutations[i,], sep = "", collapse = ">"))
      }
      i <- i + 1
    }  
    j <- j + 1
  }
  
  
  genes_of_mutations$start <- genes_of_mutations$start - range
  genes_of_mutations$end <- genes_of_mutations$end + range
  # If there is still any mutation lacking of plot it is because it is not part of the listed genes
  genes_of_mutations <- data.frame("Interval" = paste0(genes_of_mutations$chromosome,":", genes_of_mutations$start, "-", genes_of_mutations$end), "gene" = genes_of_mutations$name, "mut" = genes_of_mutations$mut, row.names = NULL)
  
  return(genes_of_mutations)
}

#' @title Designs the gene pileup plot
#' @description Checking the allelic coverage for the gene and the mutations belonging to the gene, this function designs the plot, its lines, text and legends
#' @param allelic_coverage DataFrame of colomns POS-REFCounts-ALT_counts-REF-ALT-TotalCounts
#' @param mutations_at_genes DataFrame of colomns: Interval, gene name, mutations found (REF-POS-ALT)
#' @return Recorded gene pileup plot
#' @import showtext
plot_gene_design <- function(allelic_coverage, mutations_at_genes) {
  
  Position <- allelic_coverage$POS
  Counts <- allelic_coverage$`Total Counts`
  ALT_counts <- allelic_coverage$`ALT Counts`
  par(xpd = TRUE)
  
  print(max(Counts))
  
  plot(Position, Counts, type = "l", col = "black", ylim = c(0, max(Counts)), xlab = "", ylab = "", cex.axis = 3)
  lines(Position, ALT_counts, type = "o", col = "red") # Same values for X axis
  
  title(main = mutations_at_genes$gene,
        xlab = "Position",
        ylab = "Counts",
        cex.lab = 3, cex.main = 3)
  
  # Legend
  legend("bottomleft", legend = c("Total Counts", "ALT Counts"), col = c("black", "red"), lty = 1, bty = "n", cex = 2)
  
  # Saves the mutations found in the analyzed gene
  mutations_in_gen <- mutations_at_genes$mut %>%
    strsplit(" ")
  
  # Mutations
  for (i in (1:length(mutations_in_gen[[1]]))) {
    mtext(mutations_in_gen[[1]][i], side = 3, adj = 1,
          line = (length(mutations_in_gen[[1]]) - i), cex = 2.5)
  }
  
  plot_of_gene <- recordPlot()
  
  return(plot_of_gene)
}

#' @title Designs the mutation pileup plot
#' @description Checking the allelic coverage for the mutated position and the <range> closests nucleotides, his function designs the plot, its lines, text and legends
#' @param data_mutation_plot DataFrame of colomns POS-REFCounts-ALT_counts-REF-ALT-TotalCounts for the mutated position and the <range> closests nucleotides
#' @param alternative_base ALT of the mutation
#' @return Recorded mutation pileup plot
plot_mutation_design <- function(data_mutation_plot, alternative_base, range) {
  
  # X axis: positions
  Position <- data_mutation_plot$POS
  
  # Y axis: Counts
  Counts <- data_mutation_plot$`REF Counts`
  ALT_counts <- data_mutation_plot$`ALT Counts`
  
  # ALT sequence. Only the mutation has its letter
  ALT_letters <- c(rep(".", range), alternative_base, rep(".", range))
  
  # Title: mutation
  mutation <- sprintf("%s>%s>%s", data_mutation_plot$REF[(nrow(data_mutation_plot)/2)+0.5], data_mutation_plot$POS[(nrow(data_mutation_plot)/2)+0.5], alternative_base)
  
  # First line of plot: REF sequence
  plot(Position, Counts, type = "b", pch = data_mutation_plot$REF, col = "black", cex = 3, ylim = c(-500, max(Counts)*1.1), xlab = "", ylab = "", cex.axis = 3)
  
  title(main = sprintf("%s", mutation),
        xlab = "Position",
        ylab = "Counts",
        cex.lab = 3, cex.main = 3)
  
  # Second line of plot: ALT sequence. Only the mutation has its letter
  lines(Position, ALT_counts, type = "b", pch = ALT_letters, col = "red", cex = 3)
  
  # Legend always on the right because the mutation will always be at the centre
  legend("right", legend = c("REF Counts", "ALT Counts"), col = c("black", "red"), lty = 1, bty = "n", cex = 2)
  
  plot_of_mutation <- recordPlot()
  
  return(plot_of_mutation)
}

#' @title Designs the WB that will be displayed on a XLSX file afterwards
#' @description Depending on the size of each plot, this function arranges the plot in the different cells of the WB
#' @param plots_to_xlsx Nested list of pileup plots. The main list is named after the mutated genes. Inside of them there is the gene plot and all mutation plots belonging to that gene.
#' @return WB with all the arranged plots in it
xlsx_plot_design <- function(plots_to_xlsx) {
  wb <- openxlsx::createWorkbook()
  
  # Sheets of the excel
  openxlsx::addWorksheet(wb, "Plot_Mut")
  
  #Add Plot
  for (i in 1:length(plots_to_xlsx)) { # i is related to the X axis
    for (j in 1:length(plots_to_xlsx[[i]])) { # j s related to the Y axis
      print(plots_to_xlsx[[i]][j])
      openxlsx::insertPlot(wb, 1, xy = c((i*9)-8, (j*23)-21), width = 18.5, height = 12, fileType = "png", units = "cm")
    }
  }
  
  return(wb)
}
