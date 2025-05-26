# BAM PARA IGV ####################################

archivo_bam="/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_Aligned_out.bam"
sorted_bam="/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_sorted.bam"

samtools="/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Samtools/samtools-1.16.1/samtools"
export PATH="/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Samtools/samtools-1.16.1/samtools:$PATH"
source ~/.bashrc

samtools sort -@ 8 -o "$sorted_bam" "$archivo_bam"
#indexo:
samtools index -@ 8 "$sorted_bam"

# Obtener un fragmento del bam que incluya los pbreakpoints de las fusiones que me interesan: ----
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

#-------------------------------------------------------------------------

cd /media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed

database="/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/database"
AnnotationHG38="/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.110.gtf"

/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/Arriba/arriba_v2.4.0/draw_fusions.R \
  --fusions=/media/4tb2/Daniela/Environ/FusionesActualizado/Tiroides/CarcinomasOncociticos-PRJNA445446/SRR6888826/trimmed/SRR6888826_fusions.tsv \
  --alignments=$sorted_bam \
  --output=fusions.pdf \
  --annotation=/home/juan/R/x86_64-pc-linux-gnu-library/4.1/OMICsdoSof/HG38/Homo_sapiens.GRCh38.110.gtf \
  --cytobands=$database/cytobands_hg38_GRCh38_v2.4.0.tsv \
  --proteinDomains=$database/protein_domains_hg38_GRCh38_v2.4.0.gff3
