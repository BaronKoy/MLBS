# Processes 1 pileup file per job. Requires a lit of the names of input pileup files stored in a file labelled pileup_list.txt

# import modules
import sys
import os
from collections import defaultdict
from intervaltree import IntervalTree

def load_regions(region_file):
    region_trees = defaultdict(IntervalTree)
    with open(region_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom, start, end = parts[:3]
                region_trees[chrom].addi(int(start), int(end))
    return region_trees

def match_positions(start_file, region_trees, out_file):
    with open(start_file, 'r') as infile, open(out_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom, pos = parts[0], int(parts[1])
                if chrom in region_trees and region_trees[chrom].overlaps(pos):
                    outfile.write(line)

def main():
    region_file = "/data/home/bty565/empirical/inputs/introns.txt"
    pileup_file = sys.argv[1]
    out_file = f"{os.path.splitext(pileup_file)[0]}_matches.txt"

    region_trees = load_regions(region_file)
    match_positions(pileup_file, region_trees, out_file)

if __name__ == "__main__":
    main()
