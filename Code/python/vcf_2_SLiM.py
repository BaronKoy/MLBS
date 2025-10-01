# Grabs a VCF with multiple chromosome entries and transforms into a "continuous position" file
import pandas as pd

# Define chromosome lengths
chrom_lengths = {
    "X": 23542271,
    "2L": 23513712,
    "2R": 25286936,
    "3L": 28110227,
    "3R": 32079331,
    "4": 1348131
}

# Compute cumulative offsets
chrom_order = ["X","2L","2R","3L","3R","4"]
chrom_offsets = {}
offset = 0
for chrom in chrom_order:
    chrom_offsets[chrom] = offset
    offset += chrom_lengths[chrom]

# Function to convert a VCF file
def convert_vcf(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                # Copy header lines as-is
                outfile.write(line)
            else:
                cols = line.strip().split('\t')
                chrom = cols[0]
                pos = int(cols[1])
                
                if chrom not in chrom_offsets:
                    raise ValueError(f"Chromosome {chrom} not in offsets dictionary")
                
                # Add chromosome offset to position
                new_pos = pos + chrom_offsets[chrom]
                cols[1] = str(new_pos)
                outfile.write('\t'.join(cols) + '\n')

# File name
convert_vcf('/home/baron/Documents/PhD/Data/reuter_data/VCF/dgrp2_FB.vcf', 'continuous_genome.vcf')
