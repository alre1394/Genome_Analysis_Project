#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_mapping_DNA_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/mapping_DNA.out"

#IMPORTING MODULES
module load bioinfo-tools
module load bwa/0.7.18
module load samtools/1.20

#PATHS USED
export TRIMMED_DNA="/home/alre1394/GA_PROJECT_PAPER2/data/processing/trimmed/DNA"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/genomics/assemblies"

#CONSTRUCTING FM-INDEX FOR THE REFERENCE GENOME OF STRAIN R7
bwa index "$OUTPUT/R7/assembly/assembly.fasta"

#GENOME MAPPING WITH TRIMMED DNA PAIRED-END READS (FORWARD AND REVERSE) FROM STRAIN R7
bwa mem -t 8 "$OUTPUT/R7/assembly/assembly.fasta" "$TRIMMED_DNA/SRR24413071_forward_paired.fq.gz" "$TRIMMED_DNA/SRR24413071_reverse_paired.fq.gz" | samtools view -Sb - | samtools sort -o "$OUTPUT/R7/mapping/R7_sorted.bam"
samtools index "$OUTPUT/R7/mapping/R7_sorted.bam"

#CONSTRUCTING FM-INDEX FOR THE REFERENCE GENOME OF STRAIN HP126
bwa index "$OUTPUT/HP126/assembly/assembly.fasta"

#GENOME MAPPING WITH TRIMMED DNA PAIRED-END READS (FORWARD AND REVERSE) FROM STRAIN HP126
bwa mem -t 8 "$OUTPUT/HP126/assembly/assembly.fasta" "$TRIMMED_DNA/SRR24413065_forward_paired.fq.gz" "$TRIMMED_DNA/SRR24413065_reverse_paired.fq.gz" | samtools view -Sb - | samtools sort -o "$OUTPUT/HP126/mapping/HP126_sorted.bam"
samtools index "$OUTPUT/HP126/mapping/HP126_sorted.bam"

#CONSTRUCTING FM-INDEX FOR THE REFERENCE GENOME OF STRAIN DV3
bwa index "$OUTPUT/DV3/assembly/assembly.fasta"

#GENOME MAPPING WITH TRIMMED DNA PAIRED-END READS (FORWARD AND REVERSE) FROM STRAIN DV3
bwa mem -t 8 "$OUTPUT/DV3/assembly/assembly.fasta" "$TRIMMED_DNA/SRR24413080_forward_paired.fq.gz" "$TRIMMED_DNA/SRR24413080_reverse_paired.fq.gz" | samtools view -Sb - | samtools sort -o "$OUTPUT/DV3/mapping/DV3_sorted.bam"
samtools index "$OUTPUT/DV3/mapping/DV3_sorted.bam"
