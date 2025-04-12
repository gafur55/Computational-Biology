def rle(sequence):
    if not sequence:
        return []
    
    rle_result = []
    count = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i-1]:
            count += 1
        else:
            rle_result.append((count, sequence[i-1]))
            count = 1
    rle_result.append((count, sequence[-1]))  # Append the last sequence
    return rle_result

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Concatenate the lines skipping the header line (lines starting with ">")
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def main(file_path):
    sequence = read_fasta(file_path)
    rle_result = rle(sequence)
    
    # Output the RLE and the length of both the original and the RLE sequences
    print("RLE Encoding of the Sequence:")
    print(rle_result)
    print("Length of the original sequence:", len(sequence))
    print("Length of the RLE representation:", len(rle_result))

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        main(sys.argv[1])  # Read from the specified file
    else:
        print("Please provide a path to the FASTA file.")
