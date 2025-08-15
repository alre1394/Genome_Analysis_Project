#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_unflipped_comparisons_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/unflipped_comparisons.out"

#IMPORTING MODULES
module load bioinfo-tools
module load blast/2.15.0+

#PATHS USED
export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"
export OUTPUT="$ASSEMBLY/comparisons"

#GENOME ALIGNMENT OF THE POLISHED GENOMES OF STRAINS R7 AND DV3 AGAINST STRAIN HP126
blastn -query "$ASSEMBLY/R7/polished/R7_genome_polished.fasta" \
       -subject "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
       -out "$OUTPUT/R7_HP126.txt" \
       -outfmt 6

blastn -query "$ASSEMBLY/DV3/polished/DV3_genome_polished.fasta" \
       -subject "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
       -out "$OUTPUT/DV3_HP126.txt" \
       -outfmt 6
