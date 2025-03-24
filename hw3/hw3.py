from collections import defaultdict, Counter



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

def reconstruct_from_kmer(dna):
    return(string_by_genome_path(find_eulerian_path(debruijn_from_kmers(dna))))


def kmer_composition(text, k):
    composed_kmers = []
    for i in range(len(text)- k + 1):
        composed_kmers.append(text[i:i+k])
    return(sorted(composed_kmers))


def generate_k_d_mers(genome, k, d, output_file="/home/gafur/Documents/computational_biology/hw3/k_d_mers.txt"):
    k_d_mers = []
    
    # Generate all (k,d)-mers
    for i in range(len(genome) - 2 * k - d + 1):
        p1 = genome[i:i+k]  # First k-mer
        p2 = genome[i+k+d:i+2*k+d]  # Second k-mer
        pair = p1+ '|'+ p2

        k_d_mers.append(pair)

    # Write to file
    with open(output_file, "w") as f:
        # Write first line with given numbers
        f.write(f"{k} {d}\n")
        
        # Write paired k-mers in the required format
        for pair in k_d_mers:
            f.write(f"{pair}\n")
    
    print(f"(k,d)-mer composition Results saved to {output_file}")


def read_fasta(file_path):
    """Reads a genome sequence from a FASTA file and returns it as a single string."""
    sequence = []
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            if not line.startswith(">"):  # Skip FASTA header lines
                sequence.append(line.strip())  # Remove newlines and spaces
    return "".join(sequence) 


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


def find_eulerian_path(graph):
    start, end = add_imaginary_edge(graph)
    cycle = eulerian_cycle(graph, start=end)[:-1]
    for i in range(len(cycle)):
        if cycle[i] == start and cycle[(i+1) % len(cycle)] == end:
            path = cycle[i+1:] + cycle[:i+1]

    return path


def string_by_genome_path(patterns):
    string = patterns[0]
    for i in range(1, len(patterns)):
        string += patterns[i][len(patterns[i])-1]
    return string


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






def main():
    file_path = "/home/gafur/Documents/computational_biology/hw3/sequence.fasta"
    genome_sequence = read_fasta(file_path)

    dna = kmer_composition(genome_sequence, 1000)                                                                             
    reconstructed = reconstruct_from_kmer(dna)
    if reconstructed == genome_sequence:
        print("correct")
    else:
        print("incorrect")

    # k_values = [10, 100, 1000]
    # correctly_constructed_k_d = []


    # for k in k_values:
    #     for d in range(100):
    #         generate_k_d_mers(genome_sequence, k, d)

    #         #paired debruijn graph construction
    #         file = "/home/gafur/Documents/computational_biology/hw3/k_d_mers.txt"
    #         ints, *pairs = open(file).read().splitlines()       
    #         k, d = map(int, ints.split())
    #         pairs = [x.split("|") for x in pairs]
    #         paired_debruijn_graph_construction_result = string_from_paired_composition(pairs, k, d)
    #         print("Testing with k=", k, "d=", d)
    #         if paired_debruijn_graph_construction_result == genome_sequence:
    #             print("Correctly constructed sequence")
    #             correct = 'k='+str(k) + '|'+'d='+str(d)
    #             print(correct)
    #             correctly_constructed_k_d.append(correct)
    #             break
    #         else:
    #             print("Incorrectly constructed sequence")

    # print("Final Results")
    # for i in correctly_constructed_k_d:
    #     print(i)

    """I tested all 3 different combinations (k = 10, 100, 1000) of (k,d)-mers,
    while keeping d=100 constant. With k=10, sequence was not constructed correctly, 
    however, k=100,1000 were correctly constructed."""



    



if __name__ == "__main__":
    main()