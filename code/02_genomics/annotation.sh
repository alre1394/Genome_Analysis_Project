#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_annotation_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/annotation.out"

#IMPORTING MODULES
module load bioinfo-tools
module load prokka/1.45-5b58020

#PATHS USED
export PROTEIN_REF="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics"
export ASSEMBLY="$PROTEIN_REF/assemblies"

#ANNOTATION FOR FLIPPED GENOMES OF STRAINS R7 AND DV3
for strain in R7 DV3
do
prokka --outdir "$ASSEMBLY/${strain}/annotation/flipped" \
  --prefix "${strain}_annotated" \
  --proteins "$PROTEIN_REF/ATCC10970_protein_ref.faa" \
  --kingdom Bacteria \
  --genus Streptomyces \
  --species rimosus \
  --gram pos \
  "$ASSEMBLY/${strain}/polished/${strain}_genome_polished_flipped.fasta"
done

#ANNOTATION FOR UNFLIPPED GENOMES OF STRAINS R7 AND DV3
for strain in R7 DV3
do
prokka --outdir "$ASSEMBLY/${strain}/annotation/unflipped" \
  --prefix "${strain}_annotated" \
  --proteins "$PROTEIN_REF/ATCC10970_protein_ref.faa" \
  --kingdom Bacteria \
  --genus Streptomyces \
  --species rimosus \
  --gram pos \
  "$ASSEMBLY/${strain}/polished/${strain}_genome_polished.fasta"
done

#ANNOTATION FOR UNFLIPPED GENOME OF STRAIN HP126
prokka --outdir "$ASSEMBLY/HP126/annotation" \
  --prefix HP126_annotated \
  --proteins "$PROTEIN_REF/ATCC10970_protein_ref.faa" \
  --kingdom Bacteria \
  --genus Streptomyces \
  --species rimosus \
  --gram pos \
  "$ASSEMBLY/HP126/polished/HP126_genome_polished.fasta"
