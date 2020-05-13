import sys
from collections import defaultdict


def count_ref_bases(ref):

    total = 0
    with open(ref) as f:
        for line in f:
            if line[0] != ">":
                total += len(line.strip())

    return total

def read_blocks(input_file, blocks):

    with open(input_file) as f:

        for i in range(5):
            f.next()

        for line in f:
            line_split = line.strip().split()

            read_name = line_split[9]
            block_sizes = line_split[18]
            block_starts = line_split[20]

            for start, size in zip(
                block_starts.split(',')[:-1],
                block_sizes.split(',')[:-1],
            ):
                blocks[read_name].add((
                    int(start) + 1,
                    int(start) + int(size),
                ))

        return blocks


def count_aligned_bases(block_dict, seq_length):

    def check_overlap(region_1, region_2):
        pass

    total = 0

    for blocks in block_dict.values():

        # Account for padded sequence
        _blocks = []
        for block in blocks:
            if block[0] <= seq_length and block[1] <= seq_length:
                _blocks.append(block)
            if block[0] <= seq_length and block[1] > seq_length:
                _blocks.append([block[0], seq_length])
                _blocks.append([1, block[1] - seq_length])
            if block[0] > seq_length and block[1] > seq_length:
                _blocks.append([block[0] - seq_length, block[1] - seq_length])
        blocks = _blocks

        overlap = set()
        for block in blocks:
            overlap = overlap | set(range(block[0], block[1] + 1))

        total += len(overlap)

    return total

def get_aligned_bases(input_prefix, ref):

    ref_length = count_ref_bases(ref)  

    blocks = defaultdict(set)
    blocks = read_blocks(input_prefix + '.read_1.psl', blocks)
    blocks = read_blocks(input_prefix + '.read_2.psl', blocks)

    count = count_aligned_bases(blocks, ref_length)

    with open(input_prefix + '.aligned_bases.txt', 'w') as OUTPUT:
        OUTPUT.write('{}\t{}\n'.format(input_prefix,count))

if __name__ == '__main__':

    if len(sys.argv) < 3:
        sys.stdout.write("Usage: " + sys.argv[0] + "\n           <PSL file prefix>\n           <Reference sequence>\n")
        exit()

    get_aligned_bases(sys.argv[1], sys.argv[2])
