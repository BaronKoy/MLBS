#!/usr/bin/env python3
"""
make_p3_vcf_from_slim.py

Usage:
  python make_p3_vcf_from_slim.py <input.vcf> <slim_mutations.txt> <slim_genomes.txt> <out_p3.vcf>

- slim_mutations.txt: lines like:
  #OUT: 5 5 T p3 424 m1 10056 0 0.5 p-1 1 3 A
- slim_genomes.txt: lines like:
  p3:i1 41 42 43 0 ...
  p3:i1 104 105 68 ...
  (two lines per individual = two haplotypes; script groups by sample name before space)
"""
import sys
from collections import defaultdict, OrderedDict

if len(sys.argv) != 5:
    sys.exit("Usage: python make_p3_vcf_from_slim.py <input.vcf> <slim_mutations.txt> <slim_genomes.txt> <out_p3.vcf>")

IN_VCF = sys.argv[1]
SLIM_MUTS = sys.argv[2]
SLIM_GENOS = sys.argv[3]
OUT_VCF = sys.argv[4]

# ---------------------------
# 1) Parse SLiM mutation file
# ---------------------------
# We need:
#   mutID -> (pos, derived)
#   pos -> ordered list of (mutID, derived) (to build multi-ALT)
mutid_to_info = {}          # mutID -> (pos, derived)
pos_to_mutids = OrderedDict()  # pos -> list of mutIDs (ordered by appearance)
with open(SLIM_MUTS) as fh:
    for line in fh:
        line = line.strip()
        if not line or not line.startswith("#OUT:"):
            continue
        parts = line.split()
        # Example: #OUT: 5 5 T p3 424 m1 10056 0 0.5 p-1 1 3 A
        try:
            pop = parts[4]
            if pop != "p3":
                continue
            mutid = parts[6]          # e.g. m1 (but in some outputs mut ids may be numeric; keep as str)
            pos = int(parts[7])       # genomic position
            derived = parts[-1]       # last token = derived base (A/T/C/G or multi char)
        except Exception:
            continue
        mutid_to_info[mutid] = (pos, derived)
        if pos not in pos_to_mutids:
            pos_to_mutids[pos] = []
        # Keep each mutID only once per pos
        if mutid not in pos_to_mutids[pos]:
            pos_to_mutids[pos].append(mutid)

# ---------------------------
# 2) Parse SLiM genome lists -> per-sample haplotypes
# ---------------------------
# We'll gather haplotypes as list per sample: haplotypes[sample] = [hap1_mutIDs_set, hap2_mutIDs_set, ...]
# Expect two haplotypes per individual (hap1, hap2). If only one present, treat as haploid second empty.
haplotypes = defaultdict(list)  # sample -> list of sets of mutIDs
with open(SLIM_GENOS) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        sample = parts[0]   # e.g. p3:i1
        muts = parts[1:]
        # Some entries may be numeric IDs; keep as strings to match mutIDs keys
        haplotypes[sample].append(set(muts))

# Ensure every sample has exactly 2 haplotypes (if only 1, append empty)
samples = sorted(haplotypes.keys())
for s in samples:
    if len(haplotypes[s]) == 0:
        haplotypes[s] = [set(), set()]
    elif len(haplotypes[s]) == 1:
        haplotypes[s].append(set())
    elif len(haplotypes[s]) > 2:
        # If more than 2 haplotypes present, keep first two and warn
        haplotypes[s] = haplotypes[s][:2]

# ---------------------------
# 3) Build quick mapping: pos -> derived_alleles list (in the order we will place ALTs)
# ---------------------------
pos_to_alleles = OrderedDict()   # pos -> list of derived bases (ALT order)
pos_mutid_to_allele_index = {}   # (pos, mutid) -> allele index (0-based relative to ALT list)
for pos, mutids in pos_to_mutids.items():
    alleles = []
    for mutid in mutids:
        derived = mutid_to_info.get(mutid, (None, None))[1]
        if derived is None:
            continue
        if derived not in alleles:
            alleles.append(derived)
    pos_to_alleles[pos] = alleles
    for mid in mutids:
        # assign index if derived present
        derived = mutid_to_info[mid][1]
        if derived in alleles:
            pos_mutid_to_allele_index[(pos, mid)] = alleles.index(derived)

# ---------------------------
# 4) For every sample, compute genotype per pos:
#    sample_geno[sample][pos] = (allele_index_count_on_hap1, allele_index_count_on_hap2)
#    We'll represent genotype as tuple of two allele indices: 0 => REF, 1=>first ALT, 2=>second ALT, ...
# ---------------------------
sample_genotype = {s: {} for s in samples}

for s in samples:
    hap1, hap2 = haplotypes[s][0], haplotypes[s][1]
    # For efficiency, create set of mutIDs per hap
    for pos, alleles in pos_to_alleles.items():
        # defaults: both ref
        a1 = 0
        a2 = 0
        # check hap1
        for mutid in pos_to_mutids.get(pos, []):
            if mutid in hap1:
                idx = pos_mutid_to_allele_index.get((pos, mutid), None)
                if idx is not None:
                    a1 = idx + 1
                    break
        # check hap2
        for mutid in pos_to_mutids.get(pos, []):
            if mutid in hap2:
                idx = pos_mutid_to_allele_index.get((pos, mutid), None)
                if idx is not None:
                    a2 = idx + 1
                    break
        sample_genotype[s][pos] = (a1, a2)

# ---------------------------
# 5) Write output VCF: read input header, write header but include only p3 samples (our samples list)
#    and compute AF per pos from sample_genotype
# ---------------------------
def write_vcf():
    with open(IN_VCF) as vin, open(OUT_VCF, "w") as vout:
        for line in vin:
            if line.startswith("##"):
                vout.write(line)
                continue
            if line.startswith("#CHROM"):
                # write new header with only our sample names
                head = line.strip().split("\t")[:9]
                vout.write("\t".join(head + samples) + "\n")
                continue
            # data lines
            parts = line.strip().split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            id_field = parts[2]
            ref = parts[3]
            qual = parts[5] if len(parts) > 5 else "."
            filt = parts[6] if len(parts) > 6 else "."
            info = parts[7] if len(parts) > 7 else "."
            fmt = parts[8] if len(parts) > 8 else "GT"

            if pos in pos_to_alleles:
                alts = pos_to_alleles[pos]
                alt_field = ",".join(alts) if alts else "."
                # compute AF from sample_genotype
                alt_counts = [0] * len(alts)
                total_alleles = 2 * len(samples)
                sample_gt_strings = []
                for s in samples:
                    a1, a2 = sample_genotype[s].get(pos, (0,0))
                    # count alt alleles
                    if a1 > 0:
                        alt_counts[a1-1] += 1
                    if a2 > 0:
                        alt_counts[a2-1] += 1
                    # compose GT string
                    sample_gt_strings.append(f"{a1}/{a2}")
                # AF per ALT
                afs = [round(c/total_alleles,6) for c in alt_counts]
                info_field = f"AF={','.join(map(str,afs))}"
                # prepare full line
                # minimal other fields preserved from input
                out_cols = [chrom, str(pos), id_field, ref, alt_field, qual, filt, info_field, fmt] + sample_gt_strings
                vout.write("\t".join(out_cols) + "\n")
            else:
                # invariant site: keep ALT ".", AF=0.0, genotypes 0/0
                alt_field = "."
                info_field = "AF=0.0"
                sample_gt_strings = ["0/0"] * len(samples)
                out_cols = [chrom, str(pos), id_field, ref, alt_field, qual, filt, info_field, fmt] + sample_gt_strings
                vout.write("\t".join(out_cols) + "\n")

    print("Wrote p3-only VCF to:", OUT_VCF)

write_vcf()
