def suffix_array(text):
    seqs = dict((i, text[i:]) for i in range(len(text)))
    return sorted(seqs.keys(), key=lambda x: seqs[x])

def bwt(seq):
    return "".join(seq[i - 1] for i in suffix_array(seq))


def main(file):
    print(bwt(open(file).read().rstrip()))


if __name__ == "__main__":
    import sys
    main(sys.argv[1]) 