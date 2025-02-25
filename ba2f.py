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
						

dna = []                    
with open("rosalind_ba2f.txt", "r") as file:
    for line in file:
        words = line.strip().split()  # Removes '\n' and splits by space
        dna.extend(words)  # Adds all words to the list


best_result = randomized_motif_search(dna, 15, 20)
for m in range(1000):
    result = randomized_motif_search(dna, 15, 20)
    if score(result) < score(best_result):
        best_result = result


for i in best_result:
    print(i)