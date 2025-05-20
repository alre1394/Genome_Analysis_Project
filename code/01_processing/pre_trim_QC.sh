#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 01:00:00
#SBATCH -J GA_pre_trim_QC_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/pre.trim.qc.out"

#IMPORTING MODULES
module load bioinfo-tools
module load FastQC/0.11.9

#PATHS USED
export RAW="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/01_processing/QC/pre_trim"

#QC OF SHORT ILLUMINA READS BEFORE TRIMMING
fastqc -o "$OUTPUT/DNA" -t 6 "$RAW/DNA_reads/short_reads/*.fastq.gz"
fastqc -o "$OUTPUT/RNA" -t 6 "$RAW/RNA_reads/*.fastq.qz"
