
def suffix_array(text):
    seqs = dict((i, text[i:]) for i in range(len(text)))
    return sorted(seqs.keys(), key=lambda x: seqs[x])


def main(file):
    print(*suffix_array(open(file).read().rstrip()), sep=", ")


if __name__ == "__main__":
    import sys
    main(sys.argv[1])  # only input provided, output defaults to trie_output.txt
