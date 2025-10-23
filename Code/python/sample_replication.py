import sys

input_vcf = "dgrp2_FB.vcf"
output_vcf = "duplicated_FB.vcf"
n_copies = 50

with open(input_vcf) as fin, open(output_vcf, "w") as fout:
    for line in fin:
        if line.startswith("##"):
            fout.write(line)
        elif line.startswith("#CHROM"):
            parts = line.strip().split("\t")
            fixed = parts[:9]
            samples = parts[9:]
            new_samples = []
            for s in samples:
                new_samples += [f"{s}_{i+1}" for i in range(n_copies)]
            fout.write("\t".join(fixed + new_samples) + "\n")
        else:
            parts = line.strip().split("\t")
            fixed = parts[:9]
            genotypes = parts[9:]
            new_gts = []
            for gt in genotypes:
                new_gts += [gt] * n_copies
            fout.write("\t".join(fixed + new_gts) + "\n")

