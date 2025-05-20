#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH -J GA_polishing_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/code_display/polishing.out"

#IMPORTING MODULES
module load bioinfo-tools
module load Pilon/1.24

#PATHS USED
export ASSEMBLY="/home/alre1394/GA_PROJECT_PAPER2/data/genomics/assemblies"

ÂPOLISHING GENOME ASSEMBLIES
for strain in R7 HP126 DV3
do
java -Xmx16G -jar $PILON_HOME/pilon.jar \
  --genome "$ASSEMBLY/${strain}/assembly/assembly.fasta" \
  --frags "$ASSEMBLY/${strain}/mapping/${strain}_sorted.bam" \
  --output ${strain}_genome_polished \
  --outdir "$ASSEMBLY/${strain}/polished" \
  --threads 8
done
