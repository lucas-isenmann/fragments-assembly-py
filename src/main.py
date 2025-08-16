from read_fastq import read_fastq
from overlap_graph import make_graph
from independent_set import find_IS2
from assemble import assemble

file = input("Enter file path: ")
# min_overlap = input("Enter minimum overlap: ")
min_overlap = 100
fragments = read_fastq(file)
graph, adj = make_graph(fragments, int(min_overlap))
indpt_set, complement = find_IS2(adj)
assembly = assemble(graph, list(indpt_set), list(complement))
print("Fragment sequence:")
for dna in assembly:
    print("", dna)