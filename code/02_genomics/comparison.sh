#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_comparison_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/comparison.out"

#IMPORTING MODULES
module load bioinfo-tools
module load blast/2.15.0+

#PATHS USED
export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"
export OUTPUT="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/comparisons"

#GENOME ALIGNMENT FOR VISUAL COMPARISON WITH ACT
blastn -query "$ASSEMBLY/R7/polished/R7_genome_polished.fasta" -subject "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
  -out "$OUTPUT/R7_HP126_comparison.txt" -outfmt 6

blastn -query "$ASSEMBLY/DV3/polished/DV3_genome_polished.fasta" -subject "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
  -out "$OUTPUT/DV3_HP126_comparison.txt" -outfmt 6
