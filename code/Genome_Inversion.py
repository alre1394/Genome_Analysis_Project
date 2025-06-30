#Importing SeqIO from Biopython
from Bio import SeqIO

record_list_R7 = []
with open("R7_genome_polished.fasta") as reverse_stranded:
    records = SeqIO.parse("R7_genome_polished.fasta", 'fasta')
    for record in records:
        record.seq = record.seq.reverse_complement()
        record_list_R7.append(record)

with open("R7_genome_polished_flipped.fasta", 'w') as output:
    SeqIO.write(record_list_R7, output, 'fasta')

record_list_DV3 = []
with open("DV3_genome_polished.fasta") as reverse_stranded:
    records = SeqIO.parse("DV3_genome_polished.fasta", 'fasta')
    for record in records:
        record.seq = record.seq.reverse_complement()
        record_list_DV3.append(record)

with open("DV3_genome_polished_flipped.fasta", 'w') as output:
    SeqIO.write(record_list_DV3, output, 'fasta')