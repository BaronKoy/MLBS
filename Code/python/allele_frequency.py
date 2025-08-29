#!/usr/bin/env python3
import subprocess
from pathlib import Path
import csv

# -------------------------
# PATHS
# -------------------------
base_dir = Path("/data3/scratch/bty565/cages/cage_2/bams")
pileup_dir = base_dir / "pileup"
freq_dir = pileup_dir / "frequency"
ngs_out_dir = freq_dir / "ngs_output"
plot_data_dir = ngs_out_dir / "plot_data"

fasta_ref = Path("/data/home/bty565/empirical/inputs/GCA_000001215.4_full.fna")
ngsPool_script = Path("/data/home/bty565/empirical/ngsJulia/ngsPool/ngsPool.jl")

# -------------------------
# UTILS
# -------------------------
def run(cmd, **kwargs):
    """Run a shell command and raise if fails."""
    print(f"[RUNNING] {' '.join(cmd)}")
    subprocess.run(cmd, check=True, **kwargs)

# Ensure directories exist
for d in [pileup_dir, freq_dir, ngs_out_dir, plot_data_dir]:
    d.mkdir(exist_ok=True)

# -------------------------
# STEP 1: BAM → pileup
# -------------------------
for bam in base_dir.glob("*.bam"):
    out_file = pileup_dir / f"{bam.name}.pileup.gz"
    run([
        "samtools", "mpileup", str(bam),
        "--max-depth", "0", "--min-BQ", "0",
        "--fasta-ref", str(fasta_ref),
        "-o", str(out_file)
    ])

# -------------------------
# STEP 2: Extract focal region
# -------------------------
for gzfile in pileup_dir.glob("*.gz"):
    out_file = freq_dir / f"{gzfile.stem}.pileup"
    with open(out_file, "w") as fout:
        run(["cat", str(gzfile)], stdout=fout)

# -------------------------
# STEP 3: Run ngsPool (Julia)
# -------------------------
for pileup in freq_dir.glob("*.pileup"):
    out_file = ngs_out_dir / f"{pileup.name}.out.gz"
    run([
        "julia", str(ngsPool_script),
        "--fin", str(pileup),
        "--fout", str(out_file),
        "--nChroms", "96"
    ])

# -------------------------
# STEP 4: Extract pos & MAF
# -------------------------
for gzfile in ngs_out_dir.glob("*.out.gz"):
    out_file = plot_data_dir / gzfile.stem
    with open(out_file, "w") as fout:
        p1 = subprocess.Popen(["zcat", str(gzfile)], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["awk", "{print $2\",\"$11}"], stdin=p1.stdout, stdout=fout)
        p1.stdout.close()
        p2.communicate()

# -------------------------
# STEP 5: Add generation numbers
# -------------------------
generations = {
    "S??0817_": 2,
    "S??0917_": 4,
    "S??1117_": 8,
    "S??0118_": 12,
    "S??0518_": 20,
    "S??0918_": 28,
    "S??0119_": 36,
    "S??0519_": 44,
    "S??1119_": 56,
}

for pattern, gen in generations.items():
    for file in plot_data_dir.glob(pattern + "*"):
        with open(file, "r+") as f:
            lines = f.readlines()
            f.seek(0)
            for line in lines:
                f.write(line.strip() + f",{gen}\n")
            f.truncate()

# -------------------------
# STEP 6: Remove non-variant lines
# -------------------------
for file in plot_data_dir.glob("*"):
    cleaned = []
    with open(file, "r") as f:
        for line in f:
            if ",0.0," in line:
                continue
            if "position" in line or "saf_MLE" in line:
                continue
            cleaned.append(line)
    with open(file, "w") as f:
        f.writelines(cleaned)

# -------------------------
# STEP 7 (no pandas version)
# -------------------------
final_csv = base_dir / "final_plot_data.csv"

with open(final_csv, "w", newline="") as fout:
    writer = csv.writer(fout)
    # write headers
    writer.writerow(["position", "minor_allele_freq", "generation", "source_file"])

    for file in plot_data_dir.glob("*"):
        with open(file, "r") as fin:
            for line in fin:
                parts = line.strip().split(",")
                if len(parts) != 3:  # skip malformed lines
                    continue
                pos, maf, gen = parts
                # skip header-like rows
                if pos.lower() == "position" or maf.lower() == "saf_mle":
                    continue
                writer.writerow([pos, maf, gen, file.name])

print(f"[DONE] Final concatenated dataset saved as {final_csv}")
