import time
from collections import defaultdict, Counter

# --- Part 1: Generate all 5-mers from a sequence ---
def get_kmers(seq, k=5):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

# --- Part 2: Linear search to count k-mers ---
def linear_search_counts(sequence, k=5):
    kmers = get_kmers(sequence, k)
    unique_kmers = set(kmers)

    start = time.time()
    counts = {kmer: sequence.count(kmer) for kmer in unique_kmers}
    end = time.time()

    print(f"[Linear Search] Total unique 5-mers: {len(unique_kmers)}")
    print(f"[Linear Search] Execution time: {end - start:.4f} seconds")
    return counts

# --- Part 3A: Construct suffix array ---
def suffix_array(text):
    return sorted(range(len(text)), key=lambda i: text[i:])

# --- Part 3B: Burrows-Wheeler Transform ---
def bwt_transform(text):
    sa = suffix_array(text)
    return ''.join(text[i - 1] if i > 0 else text[-1] for i in sa), sa

# --- Part 3C: Create first column from BWT ---
def first_column(bwt):
    return ''.join(sorted(bwt))

# --- Part 3D: Last-to-First Mapping ---
def last_to_first_mapping(bwt):
    first = first_column(bwt)
    count_dict = defaultdict(int)
    first_occurrence = {}
    
    for c in first:
        if c not in first_occurrence:
            first_occurrence[c] = first.index(c)

    last_to_first = []
    char_count = defaultdict(int)
    for c in bwt:
        last_to_first.append(first_occurrence[c] + char_count[c])
        char_count[c] += 1
    return last_to_first

# --- Part 3E: BWT-based exact pattern matching ---
def count_pattern_bwt(pattern, bwt, last_to_first):
    top, bottom = 0, len(bwt) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            sub = bwt[top:bottom+1]
            if symbol in sub:
                top_index = bwt[:top].count(symbol)
                bottom_index = bwt[:bottom+1].count(symbol)
                top = last_to_first[top + sub.index(symbol)]
                bottom = last_to_first[bottom - sub[::-1].index(symbol)]
            else:
                return 0
        else:
            return bottom - top + 1
    return 0

# --- Part 3F: Search all 5-mers using BWT ---
def bwt_search_counts(sequence, k=5):
    text = sequence + "$"
    bwt, _ = bwt_transform(text)
    ltf = last_to_first_mapping(bwt)
    
    kmers = get_kmers(sequence, k)
    unique_kmers = set(kmers)

    start = time.time()
    counts = {kmer: count_pattern_bwt(kmer, bwt, ltf) for kmer in unique_kmers}
    end = time.time()

    print(f"[BWT Search] Total unique 5-mers: {len(unique_kmers)}")
    print(f"[BWT Search] Execution time: {end - start:.4f} seconds")
    return counts

# --- Main Program ---
def main(fasta_file):
    # Read and clean sequence
    with open(fasta_file) as f:
        lines = f.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith(">"))

    print(f"Length of input sequence: {len(sequence)}")

    print("\n--- LINEAR SEARCH ---")
    linear_counts = linear_search_counts(sequence)

    print("\n--- BWT SEARCH ---")
    bwt_counts = bwt_search_counts(sequence)

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python part3_kmer_search.py <input_fasta_file>")
    else:
        main(sys.argv[1])
