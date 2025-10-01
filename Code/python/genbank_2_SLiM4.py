from Bio import SeqIO

# Input FASTA/GenBank file
fasta_path = "/home/baron/Documents/PhD/Data/reuter_data/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"

# Chromosomes we care about (name -> accession) - specific for Reuter evole and resequence experiment. Needs to be changed for species agnostic
wanted = {
    "X":  "NC_004354.4",
    "2L": "NT_033779.5",
    "2R": "NT_033778.4",
    "3L": "NT_037436.4",
    "3R": "NT_033777.3",
    "4":  "NC_004353.4"
}

# Get lengths from FASTA
seq_lengths = {}
for rec in SeqIO.parse(fasta_path, "fasta"):
    if rec.id in wanted.values():
        seq_lengths[rec.id] = len(rec.seq)

# Build ordered chromosome list
chromosomes = [(name, acc, seq_lengths[acc]) for name, acc in wanted.items()]

# ---- Build SLiM initialize() block ----
slim_lines = []
slim_lines.append("initialize() {")
slim_lines.append('    initializeMutationType("m1", 0.5, "f", 0.0);')
slim_lines.append('    initializeGenomicElementType("g1", m1, 1.0);')

offset = 0
recomb_lines = ["position\trate"]  # header for recombination map

for name, acc, length in chromosomes:
    start = offset
    end = offset + length - 1
    slim_lines.append(f"    // {name} ({acc}, length {length})")
    slim_lines.append(f"    initializeGenomicElement(g1, {start}, {end});")
    # recombination map entry: barrier at end of this chromosome
    recomb_lines.append(f"{end}\t0.0")  
    offset += length

slim_lines.append("    initializeRecombinationRate('recomb_map.txt');")
slim_lines.append("}")

# ---- Save outputs ----
with open("genome_initialize.slim", "w") as f:
    f.write("\n".join(slim_lines))

with open("recomb_map.txt", "w") as f:
    f.write("\n".join(recomb_lines))

print("Wrote genome_initialize.slim and recomb_map.txt")
