#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 05:00:00
#SBATCH -J preprocessing_alre1394
#SBATCH --mail -type=ALL
#SBATCH --mal -user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output=trimmed.DNA.short.out

RAW='/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023'
PRE_OUT_QC='/home/alre1394/GA_PROJECT_PAPER2/data/omics_output/genomics/short_read_qc/pre_trim'
TRIMMED='/home/alre1394/GA_PROJECT_PAPER2/data/omics_output/trimmed'
POST_OUT_QC='/home/alre1394/GA_PROJECT_PAPER2/data/omics_output/genomics/short_read_qc/post_trim'

module load bioinfo-tools FastQC/0.11.9
module load bioinfo-tools trimmomatic/0.39

#Pre-trimmining QC of short reads
fastqc -o $PRE_OUT_QC/DNA/ $RAW/DNA_reads/short_reads/*.fastq.gz
fastqc -o $PRE_OUT_QC/RNA/ $RAW/RNA_reads/*.fastq.gz

#Trimming short DNA and RNA reads 
trimmomatic PE -phred33\
-trimlog "$TRIMMED/trimo/DNA/trimlog.txt"\
$RAW/DNA_reads/short_reads/*_1.fastq.gz $ RAW/DNA_reads/short_reads/*_2.fastq.gz\
"$TRIMMED/trimo/DNA/forward_paired.fq.gz" "$TRIMMED/trimo/DNA/forward_unpaired.fq.gz"\
"$TRIMMED/trimo/DNA/reverse_paired.fq.gz" "$TRIMMED/trimo/DNA/reverse_unpaired.fq.gz"\
ILLUMINACLIP:$ADAPTERS:2:30:7\
MINLEN:

trimmomatic PE -phred33\
-trimlog "$TRIMMED/trimo/RNA/trimlog.txt\
$RAW/RNA_reads/*_1.fastq.gz $RAW/RNA_reads/*_2.fastq.gz\
"$TRIMMED/trimo/RNA/*_forward_paired.fq.gz" "$TRIMMED/trimo/RNA/*_forward_unpaired.fq.gz"\
"$TRIMMED/trimo/RNA/*_reverse_paired.fq.gz" "$TRIMMED/trimo/RNA/*_reverse_unpaired.fq.gz"\
ILLUMINACLIP:$ADAPTERS:2:30:7\
MINLEN:36

#Trimming long-read ONT DNA sequences
porechop -i RAW/DNA_reads/*.fastqc.gz -> TRIMMED/pore/*.fastq

#Post-trimming quality control
fastqc -o $POST_OUT_QC/DNA/*.fq.gz $TRIMMED/trimo/DNA/*
fastqc -o $POST_OUT_QC/RNA/*.fq.gz $TRIMMED/trimo/RNA/*
