#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_mapping_RNA_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/mapping_RNA.out"

#IMPORTING MODULES
module load bioinfo-tools
module load bwa/0.7.18
module load samtools/1.20

#PATHS USED
export TRIMMED_RNA="/home/alre1394/GA_PROJECT_PAPER2/data/01_processing/trimmed/RNA"
export REF="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies/R7/polished/R7_genome_polished_flipped.fasta"
export OUTPUT="/proj/uppmax2025-3-3/nobackup/work/alre"

#CONSTRUCTING FM-INDEX WITH THE POLISHED GENOME OF STRAIN R7
bwa index $REF

#MAPPING THE TRIMMED RNA PAIRED-END READS FOR EACH SAMPLE OF STRAIN R7 TO THE REFERENCE GENOME
for code in SRR24516462 SRR24516463 SRR24516464
do 
bwa mem -t 8 $REF "$TRIMMED_RNA/${code}_forward_paired.fq.gz" "$TRIMMED_RNA/${code}_reverse_paired.fq.gz" | samtools view -Sb -@ 8 | samtools sort -@ 8 -o "$OUTPUT/R7_${code}_sorted.bam"
samtools index --threads 8 "$OUTPUT/R7_${code}_sorted.bam"
done

##MAPPING THE TRIMMED RNA PAIRED-END READS FOR EACH SAMPLE OF STRAIN HP126 TO THE REFERENCE GENOME
for code in SRR24516459 SRR24516460 SRR24516461
do 
bwa mem -t 8 $REF "$TRIMMED_RNA/${code}_forward_paired.fq.gz" "$TRIMMED_RNA/${code}_reverse_paired.fq.gz" | samtools view -Sb -@ 8 | samtools sort -@ 8 -o "$OUTPUT/HP126_${code}_sorted.bam"
samtools index --threads 8 "$OUTPUT/HP126_${code}_sorted.bam"
done

##MAPPING THE TRIMMED RNA PAIRED-END READS FOR EACH SAMPLE OF STRAIN DV3 TO THE REFERENCE GENOME
for code in SRR24516456 SRR24516457 SRR24516458
do
bwa mem -t 8 $REF "$TRIMMED_RNA/${code}_forward_paired.fq.gz" "$TRIMMED_RNA/${code}_reverse_paired.fq.gz" | samtools view -Sb -@ 8 | samtools sort -@ 8 -o "$OUTPUT/DV3_${code}_sorted.bam"
samtools index --threads 8 "$OUTPUT/DV3_${code}_sorted.bam"
done
