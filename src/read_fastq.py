from math import inf

# reads fastq files and stores read fragments in the returned list
def read_fastq(filepath):
    with open(filepath, 'r') as f:
        file = f.readlines()
        fragments = []
        for i in range(1, len(file), 4):
            fragments.append(file[i].strip("\n"))

    return fragments

# reads dna from fasta file 
def read_dna_fasta(filepath):
    with open(filepath, "r") as f:
        f.readline()
        return f.readline()


# finds the minimum, maximum, and average length of fragments
def find_min_max_avg_frag_len(filepath):
    with open(filepath, 'r') as f:
        file = f.readlines()
        min = inf
        max = 0
        avg = 0
        for i in range(0, len(file), 4):
            metadata = file[i].split(";")
            length = int(metadata[1][7:][:-2])
            avg += length
            if min > length:
                min = length
            if max < length:
                max = length

        avg /= len(file)

    return (min, max, avg)

def count_frag_with_len(length):
    with open(filepath, 'r') as f:
        file = f.readlines()
        n = 0
        for i in range(0, len(file), 4):
            metadata = file[i].split(";")
            frag_len = int(metadata[1][7:][:-2])
            if (frag_len == length):
                n += 1 

    return n


if __name__ == "__main__":
    filepath = input("Path to fastq file: ")
    # fragments = read_fastq(filepath)
    # print(len(fragments))
    # print(find_min_max_avg_frag_len(filepath))
    # print(count_frag_with_len(10000))