#!/usr/bin/env python3

# Quick script to normal positions in a vcf file against a region-based fasta

# --- Inputs ---
input_vcf = "extracted.vcf"
output_vcf = "normalized.vcf"
fasta_start = 18515500
fasta_end   = 18526500
contig_name = "3R"

# --- Processing ---
with open(input_vcf) as fin, open(output_vcf, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue

        parts = line.strip().split("\t")
        chrom, pos = parts[0], int(parts[1])

        # Process only variants in this contig & range
        if chrom == contig_name and fasta_start <= pos <= fasta_end:
            new_pos = pos - fasta_start + 1
            parts[1] = str(new_pos)
            # Optionally, rename the contig to match your FASTA name
            parts[0] = contig_name
            fout.write("\t".join(parts) + "\n")
