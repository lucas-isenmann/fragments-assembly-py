from scipy.sparse import csr_array
from scipy.sparse.csgraph import maximum_bipartite_matching
from overlap_graph import find_max_overlap


# returns a dictionary where each key is a node u in the independent set and the corresponding value is a list of nodes from the IS such
# that for each node v in the list, at least one bridge exists between u and v
def find_connections(graph, indpt_set, complement):
    bridges = {}
    for u in indpt_set:
        bridges[u] = []
        for v in indpt_set:
            if u == v:
                continue
            else:
                for frag in complement:
                    if frag in graph[u] and v in graph[frag]: # if frag connects u and v
                        bridges[u].append(v)
                        break

    return bridges


# find the orders of the independent set where each two consecutive nodes have at least one bridge
def find_orders(indpt_set, bridges):
    if len(indpt_set) == 1:
        return [indpt_set]
    
    orders = []
    for i in range(len(indpt_set)):
        current = indpt_set[i]
        remaining = indpt_set[:i] + indpt_set[i + 1:]
        for p in find_orders(remaining, bridges):
            if p[0] in bridges[current]:
                orders.append([current] + p)
                # print([current] + p)

    return orders


# doesnt work if IS has only one element
# returns a list of the possible assemblies (in terms of a sequence of indices of the fragments) using maximum bipartite matching
def assemble(graph, indpt_set, complement):
    bridges = find_connections(graph, indpt_set, complement)
    orders = find_orders(indpt_set, bridges)
    print("Orders:\n", orders, end="\n\n")
    
    n = len(orders[0])
    assembled_dna = []
    matchings = []
    for order in orders:
        # creates a bipartite graph where one partition (the rows of the matrix) contains as nodes the gaps between two consecutive nodes
        # in the independent set and the other partition (the columns) contains the nodes from the complement 
        # the edges indicate that the node from the second set forms a bridge for the gap 
        bigraph = [0] * (n - 1)
        for i in range(n - 1):
            bigraph[i] = []
            for frag in complement:
                if frag in graph[order[i]] and order[i + 1] in graph[frag]: # f is a bridge
                    bigraph[i].append(1)
                else:
                    bigraph[i].append(0)

        # print(bigraph)
        bigraph = csr_array(bigraph)
        matching = maximum_bipartite_matching(bigraph, perm_type='column')
        matchings.append(matching)

        assembled_dna.append(find_seq(order, matching, complement))

    print("Matchings:\n", matchings, end="\n\n")

    return assembled_dna

# returns the fragments in order in terms of their indices 
def find_seq(order, matching, complement):
    seq = []
    n = len(order)
    i = 0
    j = 0
    while i < n - 1:
        seq.append(order[i])
        seq.append(complement[matching[j]])
        i += 1
        j += 1

    seq.append(order[i])

    return seq

    
# returns the dna sequence as a string assembled from the fragments
def assemble_dna(fragments, seq, rc_frags, min_ovl):
    assembled_dna = ""
    n = len(seq)

    # first round requires overlap3 to decide the sense
    frag1 = fragments[seq[0]]
    frag2 = fragments[seq[1]]
    overlap1 = find_max_overlap(frag1, frag2, min_ovl)
    overlap2 = find_max_overlap(frag1, rc_frags[seq[1]], min_ovl)
    overlap3 = find_max_overlap(rc_frags[seq[0]], frag2, min_ovl)
    if overlap1 > overlap2:
        if overlap1 > overlap3:
            assembled_dna += frag1
            assembled_dna += frag2[overlap1:] 
            frag1 = frag2  
        else:
            assembled_dna += rc_frags[seq[0]]
            assembled_dna += frag2[overlap3:]
            frag1 = frag2
    else:
        if overlap2 > overlap3:
            assembled_dna += frag1
            assembled_dna += rc_frags[seq[1]][overlap2:]
            frag1 = rc_frags[seq[1]]
        else:
            assembled_dna += rc_frags[seq[0]]
            assembled_dna += frag2[overlap3:]
            frag1 = frag2

    i = 2
    while i < n:  
        frag2 = fragments[seq[i]]

        overlap1 = find_max_overlap(frag1, frag2, min_ovl)
        overlap2 = find_max_overlap(frag1, rc_frags[seq[i]], min_ovl)
        # print(seq[i - 1], seq[i], overlap1, overlap2)

        if overlap1 > overlap2:
            assembled_dna += frag2[overlap1:]
            frag1 = frag2
        else:
            assembled_dna += rc_frags[seq[i]][overlap2:]
            frag1 = rc_frags[seq[i]]

        i += 1
        
    return assembled_dna
