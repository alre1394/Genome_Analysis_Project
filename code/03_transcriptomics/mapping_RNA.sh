#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 04:00:00
#SBATCH -J GA_mapping_RNA_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/mapping_RNA.out"

module load bioinfo-tools
module load bwa/0.7.18
module load samtools/1.20

TRIMMED_RNA="/home/alre1394/GA_PROJECT_PAPER2/data/01_processing/trimmed/RNA"
ASSEMBLIES="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"
OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/03_transcriptomics/mapped_RNA"

CONSTRUCTING FM-INDEX FOR R7 FLIPPED POLISHED GENOME
bwa index "$ASSEMBLIES/R7/polished/R7_genome_polished_flipped.fasta"
#GENOME MAPPING THE TRIMMED RNA PAIRED-END READS FROM STRAIN R7 TO ITS POLISHED AND FLIPPED GENOME
for code in SRR24516462 SRR24516463 SRR24516464
do 
bwa mem -t 8 "$ASSEMBLIES/R7/polished/R7_genome_polished_flipped.fasta" "$TRIMMED_RNA/${code}_forward_paired.fq.gz" "$TRIMMED_RNA/${code}_reverse_paired.fq.gz" | samtools view -Sb -@ 8 | samtools sort -@ 8 -o "$OUTPUT/R7_flipped_${code}_sorted.bam"
samtools index --threads 8 "$OUTPUT/R7_flipped_${code}_sorted.bam"
done

#GENOME MAPPING THE TRIMMED RNA PAIRED-END READS FROM STRAIN HP126 TO THE REFERENCE GENOME OF R7
for code in SRR24516459 SRR24516460 SRR24516461
do 
bwa mem -t 8 "$ASSEMBLIES/R7/polished/R7_genome_polished_flipped.fasta" "$TRIMMED_RNA/${code}_forward_paired.fq.gz" "$TRIMMED_RNA/${code}_reverse_paired.fq.gz" | samtools view -Sb -@ 8 | samtools sort -@ 8 -o "$OUTPUT/HP126_${code}_sorted.bam"
samtools index --threads 8 "$OUTPUT/HP126_${code}_sorted.bam"
done

#GENOME MAPPING THE TRIMMED RNA PAIRED-END READS FROM STRAIN DV3 TO THE REFERENCE GENOME OF R7
for code in SRR24516456 SRR24516457 SRR24516458
do
bwa mem -t 8 "$ASSEMBLIES/R7/polished/R7_genome_polished_flipped.fasta" "$TRIMMED_RNA/${code}_forward_paired.fq.gz" "$TRIMMED_RNA/${code}_reverse_paired.fq.gz" | samtools view -Sb -@ 8 | samtools sort -@ 8 -o "$OUTPUT/DV3_flipped_${code}_sorted.bam"
samtools index --threads 8 "$OUTPUT/DV3_flipped_${code}_sorted.bam"
done
