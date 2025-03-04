from collections import defaultdict, Counter

def parce_adj_list(in_file):
    """Parse the text file to a
    formatted graph
    """
    graph = defaultdict(list)
    for line in in_file:
        u, vs = line.strip().split(' -> ')
        u, vs = int(u), list(map(int, vs.split(',')))
        graph[u].extend(vs)
    
    return graph

def find_cycle(graph, start):
    cycle = []
    u = graph[start].pop()
    while u != start:
        cycle.append(u)
        u = graph[u].pop()
    cycle.append(u)

    toRemove = [k for k, v in graph.items() if not v]
    for k in toRemove:
        del graph[k]

    return cycle


#BA3F
def eulerian_cycle(graph, start=0):
    cycle = [start] + find_cycle(graph, start)
    updated = True  
    while updated:
        updated = False
        for i, start in enumerate(cycle):
            if start in graph:
                updated = True
                cycle = cycle[:i+1] + find_cycle(graph, start) + cycle[i+1:]
                break

    return cycle

# with open('rosalind_ba3f.txt') as inFile:
#     graph = parce_adj_list(inFile)

# with open('rosalind_BA3F_out.txt', 'w') as outFile:
#     print('->'.join(map(str, eulerian_cycle(graph))), file=outFile)



#BA3G
def add_imaginary_edge(graph):
    outgoingEdgeCounts, incomingEdgeCounts = Counter(), Counter()
    for u in graph:
        outgoingEdgeCounts[u] += len(graph[u])
        for v in graph[u]:
            incomingEdgeCounts[v] += 1 

    start = list((incomingEdgeCounts - outgoingEdgeCounts).keys())[0]
    end = list((outgoingEdgeCounts - incomingEdgeCounts).keys())[0]
    # Add imaginary edge.
    if start not in graph:
        graph[start] = []
    graph[start].append(end)
    return start, end

def find_eulerian_path(graph):
    start, end = add_imaginary_edge(graph)
    cycle = eulerian_cycle(graph, start=end)[:-1]
    for i in range(len(cycle)):
        if cycle[i] == start and cycle[(i+1) % len(cycle)] == end:
            path = cycle[i+1:] + cycle[:i+1]

    return path


# with open('/home/gafur/Documents/computational_biology/Lab7/rosalind_ba3g.txt') as inFile:
#     graph = parce_adj_list(inFile)

# with open('/home/gafur/Documents/computational_biology/Lab7/rosalind_ba3g_out.txt', 'w') as outFile:
#     print('->'.join(map(str, find_eulerian_path(graph))), file=outFile)


# BA3H
def debruijn_from_kmers(patterns):
    result = {}
    for i in patterns:
        if i[:-1] not in result:
            result[i[:-1]] = [i[1:]]  
        else:
            result[i[:-1]].append(i[1:])  
    return result

def string_by_genome_path(patterns):
    string = patterns[0]
    for i in range(1, len(patterns)):
        string += patterns[i][len(patterns[i])-1]
    return string

def reconstruct_from_kmer():
    file = "/home/gafur/Documents/computational_biology/Lab7/rosalind_ba3h.txt"
    k, *dna = open(file).read().splitlines()

    return(string_by_genome_path(find_eulerian_path(debruijn_from_kmers(dna))))

# print(reconstruct_from_kmer())


#BA3J
def dbru_paired(pairs):
    g = defaultdict(list)
    for x in pairs:
        p = tuple([x[0][:-1], x[1][:-1]])
        s = tuple([x[0][1:], x[1][1:]])
        g[p].append(s)
    return g


def string_from_paired_composition(pairs, k, d):
    path = find_eulerian_path(dbru_paired(pairs))
    a = string_by_genome_path([x[0] for x in path])
    b = string_by_genome_path([x[1] for x in path])
    return a + b[-(k + d) :]


file = "/home/gafur/Documents/computational_biology/Lab7/rosalind_ba3j.txt"
ints, *pairs = open(file).read().splitlines()
k, d = map(int, ints.split())
pairs = [x.split("|") for x in pairs]
print(string_from_paired_composition(pairs, k, d))