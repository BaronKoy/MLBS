#!/usr/bin/env python3
import subprocess
from pathlib import Path
import pandas as pd

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
        run(
            ["cat", str(gzfile)],
            stdout=subprocess.PIPE
        )

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
    out_file = plot_data_dir / gzfile.stem  # remove .gz in output filename
    with open(out_file, "w") as fout:
        run(["zcat", str(gzfile)], stdout=subprocess.PIPE)
        # awk-like filter: only position + minor allele freq
        run(
            ["awk", "{print $2\",\"$11}"],
            stdin=gzfile.open("rb"),
            stdout=fout
        )

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
# STEP 7: Concatenate all into one CSV
# -------------------------
all_data = []
for file in plot_data_dir.glob("*"):
    try:
        df = pd.read_csv(
            file,
            header=None,
            names=["position", "minor_allele_freq", "generation"]
        )
        # Drop any bad header lines that slipped through
        df = df[pd.to_numeric(df["position"], errors="coerce").notnull()]
        df["source_file"] = file.name
        all_data.append(df)
    except Exception as e:
        print(f"[WARN] Skipping {file}, error: {e}")

if all_data:
    combined = pd.concat(all_data, ignore_index=True)
    final_csv = base_dir / "final_plot_data.csv"
    combined.to_csv(final_csv, index=False)
    print(f"[DONE] Final concatenated dataset saved as {final_csv}")
else:
    print("[WARN] No files found to concatenate.")
