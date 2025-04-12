import time
import sys
import os

def suffix_array(text):
    seqs = dict((i, text[i:]) for i in range(len(text)))
    return sorted(seqs.keys(), key=lambda x: seqs[x])

def bwt(seq):
    return "".join(seq[i - 1] for i in suffix_array(seq))

def main(input_file, output_file="bwt_output.txt"):
    text = open(input_file).read().rstrip()
    result = bwt(text)
    
    with open(output_file, "w") as f:
        f.write(result)
    
    print(f"BWT result saved to '{output_file}'")

if __name__ == "__main__":
    start_time = time.time()
    
    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])  # input + output
    else:
        main(sys.argv[1])  # input only
    
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.4f} seconds")
