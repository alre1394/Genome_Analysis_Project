#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_REF_comparisons_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/REF_comparisons.out"

#IMPORTING MODULES
module load bioinfo-tools
module load blast/2.15.0+

#PATHS USED
export REFERENCE="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023/reference_genome"
export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"
export OUTPUT="$ASSEMBLY/comparisons"

#GENOME ALIGNMENT OF REFERENCE GENOMES TO POLISHED ASSEMBLY FOR VISUAL COMPARISON WITH ACT AND IDENTIFICATION OF 
#POLISHED GENOMES WITH CORRECT ORIENTATIONS RELATIVE TO REFERENCE
blastn -query "$ASSEMBLY/R7/polished/R7_genome_polished.fasta" \
       -subject "$REFERENCE/R7_genome.fasta" \
       -out "$OUTPUT/REF_R7.txt" \
       -outfmt 6

blastn -query "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta" \
       -subject "$REFERENCE/HP126_genome.fasta" \
       -out "$OUTPUT/REF_HP126.txt" \
       -outfmt 6

blastn -query "$ASSEMBLY/DV3/polished/DV3_genome_polished.fasta" \
       -subject "$REFERENCE/DV3_genome.fasta" \
       -out "$OUTPUT/REF_DV3.txt" \
       -outfmt 6
