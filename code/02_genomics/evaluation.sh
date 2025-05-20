#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_evaluation_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/evaluation.out"

#IMPORTING MODULES
module load bioinfo-tools
module load quast/5.0.2

#PATHS USED
export REF_GEN="/proj/uppmax2025-3-3/Genome_Analysis/2_Beganovic_2023/reference_genome"
export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/02_genomics/assemblies"

#EVALUATION OF UNPOLISHED ASSEMBLIES
for strain in R7 HP126 DV3
do
python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py "$ASSEMBLY/${strain}/assembly/assembly.fasta" \
  --threads 8 \
  -r "${REF_GEN}/${strain}_genome.fasta" \
  -o "$ASSEMBLY/${strain}/evaluation/unpolished"
done

#EVALUATION OF POLISHED ASSEMBLIES
for strain in R7 HP126 DV3
do
python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py "$ASSEMBLY/${strain}/polished/${strain}_genome_polished.fasta" \
  --threads 8 \
  -r "${REF_GEN}/${strain}_genome.fasta" \
  -o "$ASSEMBLY/${strain}/evaluation/polished"
done
