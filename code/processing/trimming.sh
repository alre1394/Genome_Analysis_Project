#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 04:00:00
#SBATCH -J GA_trimming_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/trimo.DNA.out"

#IMPORTING MODULES
module load bioinfo-tools
module load trimmomatic/0.39

#PATHS
export RAW="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023"
export TRIMMED="/home/alre1394/GA_PROJECT_PAPER2/data/processing/trimmed"

#TRIMMING SHORT DNA ILLUMINA READS FOR STRAINS DV3 HP126 R7
for code in SRR24413065 SRR24413071 SRR24413080
do
trimmomatic PE -phred33 -threads 6\
  -trimlog "$TRIMMED/DNA/log/${code}_trimlog.txt" \
  "$RAW/DNA_reads/short_reads/${code}_1.fastq.gz" "$RAW/DNA_reads/short_reads/${code}_2.fastq.gz" \
  "$TRIMMED/DNA/${code}_forward_paired.fq.gz" "$TRIMMED/DNA/${code}_forward_unpaired.fq.gz" \
  "$TRIMMED/DNA/${code}_reverse_paired.fq.gz" "$TRIMMED/DNA/${code}_reverse_unpaired.fq.gz" \
  ILLUMINACLIP:$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa:2:30:10 \
  MINLEN:100
done

#Trimming short RNA Illumina reads for strains HP126, DV3, R7

for code in SRR24516456 SRR24516457 SRR24516458 SRR24516459 SRR24516460 SRR24516461 SRR24516462 SRR24516463 SRR24516464
do
trimmomatic PE -phred33 -threads 6 \
  -trimlog "$TRIMMED/RNA/log/${code}_trimlog.txt" \
  "$RAW/RNA_reads/${code}_1.fastq.gz" "$RAW/RNA_reads/${code}_2.fastq.gz" \
  "$TRIMMED/RNA/${code}_forward_paired.fq.gz" "$TRIMMED/RNA/${code}_forward_unpaired.fq.gz" \
  "$TRIMMED/RNA/${code}_reverse_paired.fq.gz" "$TRIMMED/RNA/${code}_reverse_unpaired.fq.gz" \
  ILLUMINACLIP:$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE.fa:2:30:10 \
  MINLEN:36
done
