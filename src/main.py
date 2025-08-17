from read_fastq import read_fastq, read_dna_fasta
from overlap_graph import make_graph, find_RC_fragments, reverse_complement
from independent_set import find_IS2
from assemble import assemble, assemble_dna


fastq_file = input("Enter the fastq file path: ")
fasta_file = input("Enter the dna fasta file path: ")
# min_overlap = input("Enter minimum overlap: ")
min_overlap = 100

dna = read_dna_fasta(fasta_file)
fragments = read_fastq(fastq_file)
graph, adj = make_graph(fragments, int(min_overlap))
rc_frags = find_RC_fragments(fragments)
indpt_set, complement = find_IS2(adj)
assembly = assemble(graph, list(indpt_set), list(complement))

for dna_seq in assembly:
    print("Fragment sequence:\n", dna_seq)
    assembled_dna = assemble_dna(fragments, dna_seq, rc_frags, min_overlap)
    check1 = dna.find(assembled_dna)
    check2 = reverse_complement(dna).find(assembled_dna)
    if check1 >= 0:
        print("Assembly successful, starts at position:", check1, "and ends at position", check1 + len(assembled_dna))
    elif check2 >= 0:
        print("Assembly successful, starts at position:", check2, "and ends at position", check2 + len(assembled_dna))
    else:
        print("Assembly unsuccessful")

    print()