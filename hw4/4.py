def overlap_alignment(s1, s2, indel=-2, match_score=1, mismatch_score=-1):
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = j * indel
        p[j, 0] = "↑"
    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = "←"

    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            match = match_score if s1[i] == s2[j] else mismatch_score
            opt = [
                m[j, i] + match,
                m[j, i + 1] + indel,
                m[j + 1, i] + indel,
            ]
            m[new] = max(opt)
            p[new] = ["↖", "↑", "←"][opt.index(max(opt))]

    # Max score in the last column
    sc = [m[j, len(s1)] for j in range(len(s2) + 1)]
    max_score = max(sc)
    j = sc.index(max_score)
    i = len(s1)

    a1, a2 = "", ""
    while i > 0 and j > 0:
        if p[j, i] == "↖":
            a1 += s1[i - 1]
            a2 += s2[j - 1]
            j -= 1
            i -= 1
        elif p[j, i] == "←":
            a1 += s1[i - 1]
            a2 += "-"
            i -= 1
        elif p[j, i] == "↑":
            a1 += "-"
            a2 += s2[j - 1]
            j -= 1

    return max_score, a1[::-1], a2[::-1]


def main():
    seq1 = "AGCTTACGGTAGCTGATGTC"
    seq2 = "TGCAGGTACCTTGGACAGTA"
    
    
    k = 10  # overlap length

    suffix1 = seq1[-k:]
    prefix2 = seq2[:k]
    suffix2 = seq2[-k:]
    prefix1 = seq1[:k]

    # Scoring schemes
    scoring_schemes = {
        "a": {"indel": -5, "mismatch": -1, "match": 1},
        "b": {"indel": -1, "mismatch": -5, "match": 1},
        "c": {"indel": -1, "mismatch": -1, "match": 5},
    }

    for label, scores in scoring_schemes.items():
        print(f"\nScoring Scheme {label}: Indel={scores['indel']}, Mismatch={scores['mismatch']}, Match={scores['match']}")
        
        # Suffix of seq1 vs prefix of seq2
        score1, a1, a2 = overlap_alignment(suffix1, prefix2,
                                        indel=scores['indel'],
                                        mismatch_score=scores['mismatch'],
                                        match_score=scores['match'])

        # Suffix of seq2 vs prefix of seq1
        score2, b1, b2 = overlap_alignment(suffix2, prefix1,
                                        indel=scores['indel'],
                                        mismatch_score=scores['mismatch'],
                                        match_score=scores['match'])

        print(f"\nSuffix of seq1 vs Prefix of seq2:")
        print(f"  Score: {score1}")
        print(f"  Alignment Length: {len(a1)}")
        print(f"  Matches: {sum(1 for x, y in zip(a1, a2) if x == y)}")
        print(f"  Mismatches: {sum(1 for x, y in zip(a1, a2) if x != y and x != '-' and y != '-')}")
        print(f"  Indels: {a1.count('-') + a2.count('-')}")

        print(f"\nSuffix of seq2 vs Prefix of seq1:")
        print(f"  Score: {score2}")
        print(f"  Alignment Length: {len(b1)}")
        print(f"  Matches: {sum(1 for x, y in zip(b1, b2) if x == y)}")
        print(f"  Mismatches: {sum(1 for x, y in zip(b1, b2) if x != y and x != '-' and y != '-')}")
        print(f"  Indels: {b1.count('-') + b2.count('-')}")
    # seq1 = str(result[1])
    # seq2 = str(result[2]) 
    
    # num_mismatch = 0
    # if len(seq1) == len(seq2):
    #     for i in range(len(seq1)):
    #         if seq1[i] != seq2[i]:
    #             num_mismatch += 1
    # else:
    #     print("sequence length are different\n")
    # print("mismatch number: ", num_mismatch)
    # print("match number: ", len(seq1) - num_mismatch)
    # print("number of indels: ", seq1.count('-') + seq2.count('-'))

    # print(result)


if __name__ == "__main__":
    main()