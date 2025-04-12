from collections import defaultdict
import time 
import sys


def edge_match(graph, node, label):
    m = [x["n"] for x in graph[node] if x["l"] == label]
    return m[0] if m else None


def trie(seqs):
    graph = defaultdict(list)
    count = 1
    for seq in seqs:
        curr_node = 0
        for s in seq:
            match = edge_match(graph, curr_node, s)
            if match:
                curr_node = match
            else:
                graph[curr_node].append({"n": count, "l": s})
                curr_node = count
                count += 1
    return graph


def main(input_file, output_file="trie_output.txt"):
    seqs = open(input_file).read().splitlines()
    g = trie(seqs)
    with open(output_file, "w") as f:
        for n1, v in g.items():
            for n2 in v:
                f.write(f"{n1}->{n2['n']}:{n2['l']}\n")


if __name__ == "__main__":
    start_time = time.time()  # Start the timer

    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        main(sys.argv[1])

    end_time = time.time()  # End the timer
    print(f"Execution time: {end_time - start_time:.4f} seconds")