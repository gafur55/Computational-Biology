from importlib import resources
import itertools


def pam250():
    path = '/home/gafur/Documents/computational_biology/hw4/blosum.txt'
    lines = open(path).read().splitlines()
    header = lines[0].split()
    return dict([x[0], dict(zip(header, map(int, x.split()[1:])))] for x in lines[1:])


def local_alignment(s1, s2, penalty=-5):
    score = pam250()
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = 0
        p[j, 0] = "↑"
    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = "←"

    m[0, 0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            opt = [
                m[j, i] + score[s1[i]][s2[j]],
                m[j, i + 1] + penalty,
                m[j + 1, i] + penalty,
                0,
            ]
            m[new] = max(opt)
            p[new] = ["↖", "↑", "←", "↖"][opt.index(max(opt))]

    max_score = max(x for x in m.values())
    j, i = [k for k, v in m.items() if v == max_score][0]
    a1, a2 = "", ""
    while i > 0 or j > 0:
        if m[j, i] == 0:
            break
        if p[j, i] == "↖":
            a1 += s1[i - 1]
            a2 += s2[j - 1]
            j, i = j - 1, i - 1
        elif p[j, i] == "←":
            a1 += s1[i - 1]
            a2 += "-"
            i = i - 1
        elif p[j, i] == "↑":
            a1 += "-"
            a2 += s2[j - 1]
            j = j - 1

    return max_score, a1[::-1], a2[::-1]


def main():
    sequences_dict = {
        'Homo sapiens':'MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE',
        'Pan troglodytes':'MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE',
        'Gallus gallus':'MGDIEKGKKIFVQKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAEGFSYTDANKNKGKCKAAFLSHWWFLVIVRSTETCRKQNNMNDQKQMYSFHFEVTACFFSCFFFYLFSCQVSLGVRIL',
        'Danio rerio':'MGDVEKGKKVFVQKCAQCHTVENGGKHKVGPNLWGLFGRKTGQAEGFSYTDANKSKGIVWGEDTLMEYLENPKKYIPGTKMIFAGIKKKGERADLIAYLKSATS',
        'Saccharomyces cerevisiae':'MTEFKAGSAKKGATLFKTRCLQCHTVEKGGPHKVGPNLHGIFGRHSGQAEGYSYTDANIKKNVLWDENNMSEYLTNPKKYIPGTKMAFGGLKKEKDRNDLITYLKKACE'
    }

    
    for (name1, s1), (name2, s2) in itertools.combinations(sequences_dict.items(), 2):
        print(f"{name1} vs {name2}")
        result = local_alignment(s1, s2)

        seq1 = str(result[1])
        seq2 = str(result[2])
    
        num_mismatch = 0
        if len(seq1) == len(seq2):
            for i in range(len(seq1)):
                if seq1[i] != seq2[i]:
                    num_mismatch += 1
        else:
            print("sequence length are different\n")
        print("length of the alignment: ", len(seq1))
        print("mismatch number: ", num_mismatch)
        print("match number: ", len(seq1) - num_mismatch)
        print("number of indels: ", seq1.count('-') + seq2.count('-'))

        print(result)

        print()


if __name__ == "__main__":
    main()