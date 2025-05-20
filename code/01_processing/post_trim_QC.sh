#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 02:00:00
#SBATCH -J GA_post_trim_QC_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/post.trim.qc.out"

#IMPORTING MODULES
module load bioinfo-tools
module load FastQC/0.11.9

#PATHS USED
export TRIMMED="/home/alre1394/GA_PROJECT_PAPER2/data/01_processing/trimmed"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/01_processing/QC/post_trim"

#QC OF SHORT ILLUMINA READS AFTER TRIMMING
fastqc -o "$OUTPUT/DNA" -t 6 "$TRIMMED/DNA/*.fq.gz"
fastqc -o "$OUTPUT/RNA" -t 6 "$TRIMMED/RNA/*.fq.gz"
