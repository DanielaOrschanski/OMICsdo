# BAM PARA IGV ####################################
# Obtener un fragmento del bam que incluya los pbreakpoints de las fusiones que me interesan:
fs <- list.files("/media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/trimmed", full.names = TRUE, recursive = FALSE)
fs <- fs[grepl("SRR6853678_sorted.bam.tmp.", fs)]
file.remove(fs)

archivo_bam="/media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/trimmed/SRR6853678_Aligned_out.bam"
sorted_bam="/media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/trimmed/SRR6853678_sorted.bam"

samtools="/home/guadalupe/R/x86_64-pc-linux-gnu-library/4.3/OMICsdoSof/Samtools"
export PATH="/home/guadalupe/R/x86_64-pc-linux-gnu-library/4.3/OMICsdoSof/Samtools/samtools-1.16.1:$PATH"
source ~/.bashrc

samtools sort -@ 8 -o "$sorted_bam" "$archivo_bam"
#indexo:
samtools index -@ 8 "$sorted_bam"

#Armo regiones al rededor de los breakpoints de cada gen de cada fusion:
samtools view -b "$sorted_bam" 8:49118244-49119244 > /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/DYM_fusion1.bam
samtools view -b "$sorted_bam" 8:49331364-49332364 > /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/DYM_fusion2.bam

# Región en chr17 (para ambas fusiones)
samtools view -b "$sorted_bam" 17:39692761-39693761 > /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/ERBB2_fusion1.bam
samtools view -b "$sorted_bam" 17:39694488-39695488 > /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/ERBB2_fusion2.bam

# Unir todos los BAM en uno solo
samtools merge /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/fusion_gastric_ERBB2.bam /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/DYM_fusion1.bam /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/DYM_fusion2.bam /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/ERBB2_fusion1.bam /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/ERBB2_fusion2.bam

# Indexar para visualizar en IGV
samtools index /media/respaldo8t/Daniela/Gastrico-PRJNA438844-6p/SRR6853678/fusion_gastric_ERBB2.bam




./draw_fusions.R \
  --fusions=/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888689/trimmed/SRR6888689_fusions.tsv \
  --alignments=Aligned.sortedByCoord.out.bam \
  --output=fusions.pdf \
  --annotation=GENCODE19.gtf \
  --cytobands=database/cytobands_hg19_hs37d5_GRCh37_v2.5.0.tsv \
  --proteinDomains=database/protein_domains_hg19_hs37d5_GRCh37_v2.5.0.gff3
