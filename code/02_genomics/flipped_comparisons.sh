#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_flipped_comparisons_alre1394
#SBATCH --mail-type=ALL
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/flipped_comparisons.out"

#IMPORTING MODULES
module load bioinfo-tools
module load blast/2.15.0+

#PATHS USED
export REFERENCE="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023/reference_genome"
export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"
export OUTPUT="$ASSEMBLY/comparisons"

#GENOME ALIGNMENT OF THE POLISHED, FLIPPED GENOMES OF STRAINS R7 AND THE POLISHED, UNFLIPPED GENOME OF STRAIN HP126
blastn -query "$ASSEMBLY/R7/polished/R7_genome_polished_flipped.fasta" \
       -subject "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
       -out "$OUTPUT/R7_flipped_HP126.txt" \
       -outfmt 6

#GENOME ALIGNMENT OF THE POLISHED, FLIPPED GENOMES OF STRAINS DV3 AND THE POLISHED, UNFLIPPED GENOME OF STRAIN HP126
blastn -query "$ASSEMBLY/DV3/polished/DV3_genome_polished_flipped.fasta" \
       -subject "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
       -out "$OUTPUT/DV3_flipped_HP126.txt" \
       -outfmt 6
