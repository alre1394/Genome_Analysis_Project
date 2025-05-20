#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_annotation_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/annotation.out"

module load bioinfo-tools
module load prokka/1.45-5b58020

export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"

#Annotation
for strain in R7 HP126 DV3
do
prokka --outdir "$ASSEMBLY/${strain}/annotation" \
  --prefix "${strain}_annotated" \
  --kingdom Bacteria \
  --genus Streptomyces \
  --species rimosus
  --cpus 4 \
  "$ASSEMBLY/${strain}/polished/${strain}_genome_polished.fasta"
done
