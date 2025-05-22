#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -J GA_assembly_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/assembly.out"

#IMPORTING MODULES
module load bioinfo-tools
module load Flye/2.9.5

#PATHS USED
export LONG_READS="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023/DNA_reads"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"

#ASSEMBLYING GENOME WITH LONG ONT SEQUENCES FROM STRAIN R7
flye --nano-raw "$LONG_READS/SRR24413072.fastq.gz" \
  --out-dir "$OUTPUT/R7/assembly" \
  --threads 8 \
  --resume

#ASSEMBLYING GENOME WITH LONG ONT SEQUENCES FROM STRAIN HP126
flye --nano-raw "$LONG_READS/SRR24413066.fastq.gz" \
  --out-dir "$OUTPUT/HP126/assembly" \
  --threads 8 \
  --resume

#ASSEMBLYING GENOME WITH LONG ONT SEQUENCES FROM STRAIN DV3
flye --nano-raw "$LONG_READS/SRR24413081.fastq.gz" \
  --out-dir "$OUTPUT/DV3/assembly" \
  --threads 8 \
  --resume
