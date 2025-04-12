import time
import sys


def suffix_array(text):
    seqs = dict((i, text[i:]) for i in range(len(text)))
    return sorted(seqs.keys(), key=lambda x: seqs[x])


def main(input_file, output_file="suffix_array_output.txt"):
    result = suffix_array(open(input_file).read().rstrip())
    with open(output_file, "w") as f:
        f.write(", ".join(map(str, result)))
    print(f"Output written to {output_file}")


if __name__ == "__main__":
    start_time = time.time()  # Start timer

    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        main(sys.argv[1])

    end_time = time.time()  # End timer
    print(f"\nExecution time: {end_time - start_time:.4f} seconds")
