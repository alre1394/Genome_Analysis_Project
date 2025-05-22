#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 06:00:00
#SBATCH -J GA_read_counting_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/read_count.out"

module load bioinfo-tools
module load htseq/2.0.2

export ASSEMBLIES="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"
export MAPPED_RNA="/home/alre1394/GA_PROJECT_PAPER2/data/03_transcriptomics/mapped_RNA"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/03_transcriptomics/read_counts"

#COUNTING MAPPED RNA READS OF GENOMICS FEATURES
#FOR R7 FLIPPED GENOME
for code in SRR24516462 SRR24516463 SRR24516464
do
htseq-count --mode union -f bam -t CDS -r pos -i ID \
            "$MAPPED_RNA/R7_flipped_${code}_sorted.bam" \
            "$ASSEMBLIES/R7/annotation/flipped/R7_annotated_del.gff" > "$OUTPUT/R7/${code}_counts.txt"
done

#FOR DV3 FLIPPED GENOME
for code in SRR24516456 SRR24516457 SRR24516458
do
htseq-count --mode union -f bam -t CDS -r pos -i ID \
            "$MAPPED_RNA/DV3_flipped_${code}_sorted.bam" \
            "$ASSEMBLIES/DV3/annotation/flipped/DV3_annotated_del.gff" > "$OUTPUT/DV3/${code}_counts.txt"
done

#FOR HP126 UNFLIPPED GENOME
for code in SRR24516459 SRR24516460 SRR24516461
do
htseq-count --mode union -f bam -t CDS -r pos -i ID\
            "$MAPPED_RNA/HP126_${code}_sorted.bam" \
            "$ASSEMBLIES/HP126/annotation/HP126_annotated_del.gff" > "$OUTPUT/HP126/${code}_counts.txt"
done
