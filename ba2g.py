import random



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

  
k = 15
t = 20
n = 2000
dna = []

with open("rosalind_ba2g.txt", "r") as file:
    for line in file:
        words = line.strip().split()  # Removes '\n' and splits by space
        dna.extend(words)  # Adds all words to the list

best = gibbs_sampler(dna, k, t, n)
s = score(best)
print(s)
for x in range(20):
    sample = gibbs_sampler(dna, k, t, n)
    print(score(sample))
    if score(sample) < s:
        s = score(sample)
        best = sample[:]

for b in best:
	print(b)

