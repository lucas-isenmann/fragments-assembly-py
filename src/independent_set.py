from math import inf


# both functions use a dictionary where the key is the node and the value is a list of the nodes adjacent to it (there 
# exists an edge betweeen them) 
# both return the independent set and the complement (remaining nodes in the graph not in the IS) to be used for assembly


# finds independent set by chosing the vertex with the minimum degree
def find_IS1(adj):
    indpt_set = set()
    complement = set()

    while adj != {}:
        v = find_min_degree_vertex(adj)
        indpt_set.add(v)
        for u in adj[v]:
            if (u not in indpt_set) and (u not in complement):
                complement.add(u)
                del adj[u]
            
        del adj[v]

    return (indpt_set, complement)

# finds independent set by looping through vertices instead
def find_IS2(adj):
    indpt_set = set()
    complement = set()

    n = len(adj)
    for v in range(n):
        if v in adj:
            indpt_set.add(v)
            for u in adj[v]:
                if (u not in indpt_set) and (u not in complement):
                    complement.add(u)
                    del adj[u]
                
            del adj[v]

    return (indpt_set, complement)


# finds minimum degree vertex using a dictionary of nodes and the nodes adjacent to them
def find_min_degree_vertex(adj):
    min = inf
    min_vertex = 0
    for v in adj.keys():
        if min > len(adj[v]):
            min = len(adj[v])
            min_vertex = v

    return min_vertex
