from read_fastq import read_fastq


# reverses fragment
def reverse_frag(frag):
    return frag[::-1]

# finds the complement of a fragment
def complement(frag):
    comp = ""
    n = len(frag)
    for i in range(n):
        if frag[i] == "A":
            comp += "T"
        elif frag[i] == "T":
            comp += "A"
        elif frag[i] == "C":
            comp += "G"
        elif frag[i] == "G":
            comp += "C"
    return comp

# returns a list of the reverse complement of each fragment 
def find_RC_fragments(fragments):
    rc_frags = []
    for frag in fragments:
        rc_frags.append(complement(reverse_frag(frag)))

    return rc_frags


# removes fragments contained within other fragments and fragments of smaller size than the min overlap
def clean_fragments(fragments, rc_frags, minOvl):
    i = 0
    n = len(fragments)
    while i < n and n != 0:
        frag1 = fragments[i]
        # first check if the fragment has a smaller length than the min overlap
        if len(frag1) < minOvl:
            fragments.remove(frag1)
            rc_frags.pop(i)
            n -= 1
            continue

        # check if any of the other fragments or their reverse complement are contained in frag1
        j = 0
        while j < n:
            frag2 = fragments[j]
            if frag1 != frag2:
                if frag2 in frag1 or rc_frags[j] in frag1 or rc_frags[j] in rc_frags[i]:
                    fragments.remove(frag2)
                    rc_frags.pop(j)
                    n -= 1
                    if j < i:
                        i -= 1
                    continue
            
            j += 1
        i += 1


#  returns overlap graph and adj, the dictionary of nodes and their neighbours (adj is later used to find independent set)
def make_graph(fragments, minOvl):
    rc_frags = find_RC_fragments(fragments)
    clean_fragments(fragments, rc_frags, minOvl)
    n = len(fragments)
    graph = {}
    adj = {}

    for i in range(n):
        graph[i] = []

        for j in range(n):
            if i == 0:
                adj[j] = set() # for each node, initialises sets in which to store the adjacent nodes
            if i != j:
                # possible cases:
                # - they overlap
                # - frag at i overlaps with reverse complement of frag at j
                # - reverse complement of frag at i overlaps with frag at j
                # - they do not overlap at all 
                overlap1 = find_max_overlap(fragments[i], fragments[j], minOvl)
                overlap2 = find_max_overlap(fragments[i], rc_frags[j], minOvl)
                overlap3 = find_max_overlap(rc_frags[i], fragments[j], minOvl)
                overlap = max(overlap1, overlap2, overlap3)
            
                if (overlap >= minOvl):
                    print("\rAdding " + str(j) + " to " + str(i), end='')
                    graph[i].append(j)
                    # the two nodes are adjacent
                    adj[i].add(j)
                    adj[j].add(i)
                    
    print()
    return (graph, adj)



# finds the max overlap between two fragments 
def find_max_overlap(frag1, frag2, minOvl):
    m = len(frag1)
    n = len(frag2)

    if (m < minOvl or n < minOvl):
        return -1
    
    # uses minimum overlap as a starting point of comparison 
    start = minOvl
    f1 = frag1[-start:]
    f2 = frag2[ :start]

    longest = 0
    if f1 == f2:
        longest = start

    for i in range(start + 1, min(m, n) + 1): 
        f1 = frag1[-i] + f1
        f2 = f2 + frag2[i - 1]    
        if f1 == f2:
            longest = i

    return longest




if __name__ == "__main__":
    fragments = read_fastq("Files/r4.fastq")
    frags = fragments.copy()
    graph, adj = make_graph(frags, 10)
    print(graph)
    # find the original indices of the fragments 
    for f in graph.keys():
        print(fragments.index(frags[f]))



'''
# make_graph test with hardcoded fragments
if __name__ == "__main__":
    fragments = ['ACGCA', 'AC', 'CGCAT', 'AT', 'TCGCG', 'CATTC', 'ATTCG']
    graph, adj = make_graph(fragments, 2)
'''


'''
# find_max_overlap test 
if __name__ == "__main__":
    s1 = input()
    s2 = input()
    print(find_max_overlap(s1, s2, 2))
'''


'''
# clean_fragments test 
if __name__ == "__main__":
    # fragments = ['ACTC', 'TGCG', 'TCG', 'ACTCTG', 'GCG', 'c', 'GT']
    # fragments = ['ACTCTG', 'ACTC']
    fragments = read_fastq("Files/r4.fastq")
    frags = fragments.copy()
    rc_frags = find_RC_fragments(fragments)
    clean_fragments(frags, rc_frags, 10)
    print(len(frags))
    for f in frags:
        print(fragments.index(f))
'''
    
