import random

def profile_most_probable_kmer(dnas, k, profile):
    kmer_dict = {}
    for dna in dnas:
        kmers_probs = []
        for i in range(len(dna)-k+1):
            kmer = dna[i:i+k]
            kmer_prob = 1
            for i in range(len(kmer)):
                if kmer[i] == 'A':
                    kmer_prob *= profile[0][i]
                elif kmer[i] == 'C':
                    kmer_prob *= profile[1][i]
                elif kmer[i] == 'G':
                    kmer_prob *= profile[2][i]
                elif kmer[i] == 'T':
                    kmer_prob *= profile[3][i]
            kmers_probs.append(kmer_prob)				
        kmer_dict[dna] = kmers_probs
    print(kmer_dict)
            
        


def n_to_s(x):
	if x == 0:
		return "A"
	if x == 1:
		return "C"
	if x == 2:
		return "G"
	if x == 3:
		return "T"


def s_to_n(symbol):
    if symbol == "A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    if symbol == "T":
        return 3

def profileForm(motifs):
	k = len(motifs[0])
	profile = [[1 for i in range(k)] for j in range(4)]
	for x in motifs:
		for i in range(len(x)):
			j = s_to_n(x[i])
			profile[j][i] += 1
	for i in range(len(profile)):
		for j in range(len(profile[i])):
			profile[i][j] = profile[i][j]/(len(motifs)+4)
			
	return(profile)



def profileProbable(text, k, profile):
	maxprob = 0
	kmer = text[0:k]
	for i in range(0,len(text) - k +1):
		prob = 1
		pattern = text[i:i+k]
		for j in range(k-1):
			l = s_to_n(pattern[j])
			
			prob *= profile[l][j]
		if prob > maxprob:
			maxprob = prob
			kmer = pattern
	return kmer

def consensus(profile):
	str = ""
	for i in range(len(profile[0])):
		max = 0
		loc = 0
		for j in range(4):
			if profile[j][i] > max:
				loc = j
				max = profile[j][i]
		str+=n_to_s(loc)
	return str

def score(motifs):
	profile = profileForm(motifs)
	cons = consensus(profile)
	score = 0
	for x in motifs:
		for i in range(len(x)):
			if cons[i] != x[i]:
				score += 1
	return score


def randomized_motif_search(dna, k, t):
	
	best_motifs = []
	motifs = []
	for seq in dna:
		random_num = random.randint(0, len(dna[0])-k)
		motifs.append(seq[random_num:random_num+k])
		
	best_motifs = motifs

	while True:
		new_motifs = []
		profile = profileForm(best_motifs)
		for i in range(len(dna)):
			new_motifs.append(profileProbable(dna[i], k, profile))
		if score(new_motifs) < score(best_motifs):
			best_motifs = new_motifs
		else:
			return best_motifs
						


def profileRandom(k, profile, text):
    probs = []
    for i in range(0,len(text) - k +1):
        prob = 1.0
        pattern = text[i:i+k]
        for j in range(k):
            l = s_to_n(pattern[j])
            prob *= profile[l][j]
        probs.append(prob)
    r = myRandom(probs)
    return r

def myRandom(dist):
    s = 0.0
    for x in dist:
        s+= x
    i = random.random()
    partial = 0.0
    for x in range(len(dist)):
        partial += dist[x]
        if partial/s >= i:
            return x

def gibbs_sampler(dna, k, t, n):
    best_motifs = []
    motifs = []
    for seq in dna:
        random_num = random.randint(0, len(dna[0])-k)
        motifs.append(seq[random_num:random_num+k])
    
    best_motifs = motifs
    for i in range(n):
        j = random.randint(0, t-1)
        profile = profileForm(motifs[:j] + motifs[j+1:])
        r = profileRandom(k, profile, dna[j])
        motifs[j] = dna[j][r:r+k]
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs

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

def approx_kmers_find(pattern, text, d):
    result = []
    for i in range(len(text)-len(pattern)+1):
        if hamming_dist(text[i:i+len(pattern)], pattern) <= d:
            result.append(i)
    return result

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

def hamming_dist(str1, str2):
    hamming_distance = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hamming_distance += 1
    return hamming_distance


def approximate_occurrences_of_pattern(pattern, text, d):
    result = []
    for i in range(len(text)-len(pattern)):
        if hamming_dist(text[i:i+len(pattern)], pattern) <= d:
            result.append(i)
    return result

def main():
    dna = []                    
    with open("/home/gafur/Documents/computational_biology/hw2/DosR.txt", "r") as file:
        for line in file:
            words = line.strip().split()  # Removes '\n' and splits by space
            dna.extend(words)  # Adds all words to the list

    

# # Brute force motif finding
#     brute_force_result = list(motif_finding(dna, 15, 3))
#     for i in brute_force_result:
#          print(i)


    mc_motifs = [
    "GGACTTCAGGCCCTA",
    "GGTCAAACGACCCTA",
    "GGACGTAAGTCCCTA",
    "GGGCTTCCAACCGTG",
    "GGCCGAACGACCCTA",
    "GGACCTTCGGCCCCA",
    "GGACTTCTGTCCCTA",
    "GGACTTTCGGCCCTG",
    "GGACTAACGGCCCTC",
    "GGACGTCCGCGACGA"]
    mc_profile = profileForm(mc_motifs)
    mc_consesus = consensus(mc_profile)

    gibbs_motifs = [
    "GGGACTTCAGGCCCT",
    "GGGTCAAACGACCCT",
    "GGGACGTAAGTCCCT",
    "CGGGCTTCCAACCGT",
    "GTGACCGACGTCCCC",
    "AGGACCTTCGGCCCC",
    "GGGACTTCTGTCCCT",
    "GGGACTTTCGGCCCT",
    "AGGACTAACGGCCCT",
    "GGGACCGAAGTCCCC"]
    gibbs_profile = profileForm(gibbs_motifs)
    gibbs_consesus = consensus(gibbs_profile)


    # Cheking the occurences of consesus motif in whole genome sequence
    with open("/home/gafur/Documents/computational_biology/hw2/sequence.fasta", "r", encoding="utf-8") as file:
        mycobacterium_tuberclosis_sequence = file.read()


    # consesus motif from Monte Carlo 
    mc_occurences = approximate_occurrences_of_pattern(mc_consesus, mycobacterium_tuberclosis_sequence, 3)
    print(mc_occurences)
    print(len(mc_occurences))

    #consesus motif from Gibbs
    gibbs_occurences = approximate_occurrences_of_pattern(gibbs_consesus, mycobacterium_tuberclosis_sequence, 3)
    print(gibbs_occurences)
    print(len(gibbs_occurences))





    # print("consesus of the motifs from MC Algo: ", mc_consesus)
    # print("consesus of the motifs from Gibbs Algo: ", gibbs_consesus)

# # Gibbs sampler results
#     k = 15
#     t = 10
#     n = 2000
#     gibbs_results = gibbs_sampler(dna, k, t, n)
#     s = score(gibbs_results)
#     print(s)
#     for x in range(20):
#         sample = gibbs_sampler(dna, k, t, n)
#         # print(score(sample))
#         if score(sample) < s:
#             s = score(sample)
#             gibbs_results = sample[:]
#     for b in gibbs_results:
#         print(b)

# # Testing with Monte Carlo pseudocounts
#     mc_results = randomized_motif_search(dna, 15, 10)
#     for m in range(10000):
#         result = randomized_motif_search(dna, 15, 10)
#         if score(result) < score(mc_results):
#             mc_results = result
#     for i in mc_results:
#         print(i)



if __name__ == "__main__":
    main()