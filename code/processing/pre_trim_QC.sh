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

#Importing modules
module load bioinfo-tools
module load FastQC/0.11.9

#Paths used
export RAW="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/preprocessing/QC/pre_trim"

#QC of short reads before trimming
fastqc -o "$OUTPUT/DNA" -t 6 "$RAW/DNA_reads/short_reads/*.fastq.gz"
fastqc -o "$OUTPUT/RNA" -t 6 "$RAW/RNA_reads/*.fastq.qz"
