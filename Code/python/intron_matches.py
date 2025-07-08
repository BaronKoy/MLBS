#Module loads
import os
import glob

def load_regions(region_file):
    regions = []
    if not os.path.exists(region_file):
        print(f"ERROR: Region file '{region_file}' not found.")
        return regions

    with open(region_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                chrom, start, end = parts[:3]
                regions.append((chrom, int(start), int(end)))
    return regions

def match_positions(start_file, regions):
    matches = []
    with open(start_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom, pos = parts[0], int(parts[1])
                for r_chrom, r_start, r_end in regions:
                    if chrom == r_chrom and r_start <= pos <= r_end:
                        matches.append(line.strip())
                        break
    return matches

def main():
    region_file = "/data/home/bty565/empirical/inputs/introns.txt"
    regions = load_regions(region_file)

    if not regions:
        print("No regions loaded. Exiting.")
        return

    start_files = glob.glob("/data3/scratch/bty565/bams/test_cage/S0*.pileup")
    if not start_files:
        print("No start file(s) found.")
        return

    for start_file in start_files:
        print(f"Processing {start_file}...")
        matched = match_positions(start_file, regions)

        if matched:
            out_file = f"{os.path.splitext(start_file)[0]}_matches.txt"
            with open(out_file, 'w') as f:
                f.write('\n'.join(matched) + '\n')
            print(f"Matches written to {out_file}")
        else:
            print(f"No matches found in {start_file}")

if __name__ == "__main__":
    main()
