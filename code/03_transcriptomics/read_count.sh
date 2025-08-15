#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -J GA_read_counting_alre1394
#SBATCH --mail-type=ALL
#SBATCH --mail-user alexander-robert.renlund.1394@student.uu.se
#SBATCH --output="/home/alre1394/GA_PROJECT_PAPER2/code/00_code_display/read_count.out"

module load bioinfo-tools
module load htseq/2.0.2

DATA="/home/alre1394/GA_PROJECT_PAPER2/data"
export REF="$DATA/02_genomics/assemblies/R7/annotation/flipped/R7_annotated_del.gff"
export MAPPED_RNA="/proj/uppmax2025-3-3/nobackup/work/alre"
export OUTPUT="$DATA/03_transcriptomics/read_counts"

#COUNTING READS FROM MAPPED RNA FEATURES OF EACH SAMPLE FROM R7
for code in SRR24516462 SRR24516463 SRR24516464
do
htseq-count -f bam -r pos -s reverse -t CDS -i ID --mode union --nonunique none \
            --secondary-alignments ignore --supplementary-alignments ignore \
            "$MAPPED_RNA/R7_${code}_sorted.bam" \
            $REF > "$OUTPUT/R7/R7_${code}_counts.txt"
done

#COUNTING READS FROM MAPPED RNA FEATURES OF EACH SAMPLE FROM HP126, with R7 ANNOTATION AS REFERENCE
for code in SRR24516459 SRR24516460 SRR24516461
do
htseq-count -f bam -r pos -s reverse -t CDS -i ID --mode union --nonunique none \
            --secondary-alignments ignore --supplementary-alignments ignore \
            "$MAPPED_RNA/HP126_${code}_sorted.bam" \
            $REF > "$OUTPUT/HP126/HP126_${code}_counts.txt"
done

#COUNTING READS FROM MAPPED RNA FEATURES OF EACH SAMPLE FROM DV3, with R7 ANNOTATION AS REFERENCE
for code in SRR24516456 SRR24516457 SRR24516458
do
htseq-count -f bam -r pos -s reverse -t CDS -i ID --mode union --nonunique none \
            --secondary-alignments ignore --supplementary-alignments ignore \
            "$MAPPED_RNA/DV3_${code}_sorted.bam" \
            $REF > "$OUTPUT/DV3/DV3_${code}_counts.txt"
done
