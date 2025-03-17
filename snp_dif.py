from Bio import SeqIO
from Bio.Seq import Seq

def count_differences(seq1, seq2):
    return sum(1 for base1, base2 in zip(seq1, seq2) if base1 != base2 and base1 != "-" and base2 != "-")

fasta_file = "./core_gene_alignment.aln"

sequences = []
for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(record.seq)

matrix = [[0] * len(sequences) for _ in range(len(sequences))]

for i in range(len(sequences)):
    for j in range(i+1, len(sequences)):
        diff_count = count_differences(sequences[i], sequences[j])
        matrix[i][j] = diff_count
        matrix[j][i] = diff_count

with open("out.txt", "w") as output_file:
    for row in matrix:
        output_file.write(" ".join(str(val) for val in row))
        output_file.write("\n")

