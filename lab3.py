#hamming distance
def hamming_dist(str1, str2):
    hamming_distance = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hamming_distance += 1
    return hamming_distance

#BA1N d-neighbors
def neighbors(pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ('A','C','G','T')
    
    neighborhood = []
    suffix_neighbors = neighbors(pattern[1:len(pattern)], d)

    for text in suffix_neighbors:
        if hamming_dist(pattern[1:len(pattern)], text) < d:
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighborhood.append(nucleotide+text)
        else:
            neighborhood.append(pattern[0]+text)
    
    return set(neighborhood)

# a = (neighbors('ACG', 1))
# for i in a:
#     print(i)



# BA2A solution
def motif_finding(dna, k, d):
    patterns = set()
    all_neighbors = set()
    first_dna = dna[0]

    for string in dna:
        for i in range(len(first_dna)-k+1):
            new_neighbors = neighbors(string[i:i+k], d)
            all_neighbors |= new_neighbors

    for neighbor in all_neighbors:
        flag = True
        for i in range(0, len(dna)):
            if len(approx_kmers_find(neighbor, dna[i], d)) == 0:
                flag = False
        if flag:
            patterns.add(neighbor)
    
    return patterns
    # result = ''
    # for i in patterns:
    #     result += i
    #     result += ' '
    # print(result)


  
# solution for find approximate kmers
def approx_kmers_find(pattern, text, d):
    result = []
    for i in range(len(text)-len(pattern)+1):
        if hamming_dist(text[i:i+len(pattern)], pattern) <= d:
            result.append(i)
    return result



# lines = []
# with open("rosalind_ba2h.txt", "r") as file:
#     for line in file:
#         words = line.strip().split()  # Removes '\n' and splits by space
#         lines.extend(words)  # Adds all words to the list

# print(lines)




# BA2H
def distance_between_pattern_and_string(pattern, dna):
    k = len(pattern)
    distance = 0
    for string in dna:
        ham_dist = 10**100
        for i in range(len(string)-k+1):
            if ham_dist > hamming_dist(pattern, string[i:i+k]):
                ham_dist = hamming_dist(pattern, string[i:i+k])
        distance += ham_dist
    return distance


# print(distance_between_pattern_and_string('AAGAGCC', lines))

#BA2B
def meadian_string(dna, k):
    distance = 10**100
    all_kmers = generate_kmers(k, ['A', 'T', 'C', 'G'])

    median = []

    for i in range(len(dna)):
        for kmer in all_kmers:
            if distance > distance_between_pattern_and_string(kmer, dna):
                distance = distance_between_pattern_and_string(kmer, dna)
                median.append(kmer)
        
    return median

def generate_kmers(k, bases):
    if k == 1:
        return bases
    small_kmers = generate_kmers(k-1, bases)
    kmers = []
    for kmer in small_kmers:
        for b in bases:
            kmers.append(kmer+b)
    return kmers

# print(meadian_string(['CCCTCCCAATCCATTGTACTGCGGCCAGCCTCCACTAGTCAC', 'ACAAAGTTTGCATCCAGTCAATCCCGGCTAGTGACCCTGACA', 'GCCAAGCAAGCCCGTCTCCGTCTAGACTCTTGGTGCCTCACC', 'ATATGACAATCCTGCTAGTCCCTGGCTCAGTATACTTGGGCT', 'TAAGCTTAGTCCATGTTCAAAAAACAATCCTACGGGTCCGAA', 'GCCAAACAAACCCGAAGCCGACGGATCCATTTCAGTCGCACG', 'GAAGAGGATAGGGCTCAGTTTTGACAAGCCATTCGGCGGGCC', 'TTTGCCTAACGACAACCCTTTCAACATTCAAGTCCTCCTCGA', 'CGTACAGTAGTAAACAAGGCGCAGGCCGTGCAATCCTACGAG', 'CAAGCCGACAGTATAAAACCCCTCACAACATCCCGATTGAGT'], 6))


#BA2C
def profile_most_probable_kmer(text, k, profile):
    kmer_dict = {}
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        kmer_probability = 1
        for i in range(len(kmer)):
            if kmer[i] == 'A':
                kmer_probability *= profile[0][i]
            elif kmer[i] == 'C':
                kmer_probability *= profile[1][i]
            elif kmer[i] == 'G':
                kmer_probability *= profile[2][i]
            elif kmer[i] == 'T':
                kmer_probability *= profile[3][i]

        kmer_dict[kmer] = kmer_probability
    
    max_val = max(kmer_dict.values())
    for key in kmer_dict:
        if kmer_dict[key] == max_val:
            print(key)

profile_most_probable_kmer('ATCAACAGCAGGGTGACTGCTTAGATTATGGATGCGCGTAATCGCGATGGATTGACCGCGCTGCGGCCTTTAGATGCTTTGTGATGGCATGGCAGGTTTATTTGAAGACACCCCATTAAGTGCCGGTACTGGGTATCTTGCTATTAACTCGATAGAATGGTAGTTCTGCGCAAGGGGGCCCTACCTGGACGGGGATGCGA', 6, [[0.212, 0.061, 0.182, 0.273, 0.182, 0.303],[0.273, 0.333, 0.242, 0.303, 0.212, 0.303,
],[0.333, 0.242, 0.333, 0.273, 0.333, 0.121
],[0.182, 0.364, 0.242, 0.152, 0.273, 0.273
]])