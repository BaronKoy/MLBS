#!/usr/bin/env bash
# Convert BAM files to sync format for PoPoolation2/Poolseq
# Assumes bams are named with S01MMYY or similar format
# Assumes intron positions are in introns.bed
# Assumes reference fasta is GCA_000001215.4_full.fna
# Baron Koylass 2025
# v1.0

# TODO...integrate into pipeline with input variables for cage number

# Moule load required software
module load samtools
module load bedtools2

# change to working directory - specific cage - TODO when integrating into pipeline this should be an input variable
cd /data3/scratch/bty565/cages/cage_1 # CHANGE AS REQUIRED
mkdir sync

# Filter bams to only include neutral SNPs
for x in *.bam ; do samtools index "$x" ; done
for x in *.bam ; do bedtools intersect -abam "$x" -b /data/home/bty565/empirical/inputs/introns.bed > "$x"_intersect.bam ; done

# Rename and copy bams to sync folder
# Collect files
files=(S01*.bam_intersect.bam) # CHANGE PREFIX AS REQUIRED

declare -A date_map

for file in "${files[@]}"; do
  # Extract last 4 digits after the initial prefix
  mm_yy=$(echo "$file" | sed -E 's/^S01([0-9]{4}).*/\1/') # CHANGE PREFIX AS REQUIRED
  
  # Parse month/year
  mm=${mm_yy:0:2}
  yy=${mm_yy:2:2}
  yyyy="20$yy"
  sortable="${yyyy}-${mm}"
  
  date_map["$file"]=$sortable
done

count=0
for file in "${!date_map[@]}"; do
  echo "${date_map[$file]} $file"
done | sort -k1,1 | while read -r date file; do
  month=$(date -d "$date-01" +%b)
  year=$(date -d "$date-01" +%Y)
  ((count++))
  
  echo "$month $year > ${count}.bam"
  cp "$file" sync/"${count}.bam"
done

# Filter bams to only include reads with MAPQ >= 1 and remove duplicates
cd sync
for x in *.bam ; do samtools view -b -q 1 -F 1024 "$x" > "$x".filtered.bam ; done
# Index filtered bams
for x in *.filtered.bam ; do samtools index "$x" ; done
# Generate mpileup
samtools mpileup 1.bam.filtered.bam 2.bam.filtered.bam 3.bam.filtered.bam 4.bam.filtered.bam 5.bam.filtered.bam 6.bam.filtered.bam 7.bam.filtered.bam 8.bam.filtered.bam 9.bam.filtered.bam --max-depth 0 --min-BQ 0 --fasta-ref /data/home/bty565/empirical/inputs/GCA_000001215.4_full.fna -o cage_1.mpileup
# Convert mpileup to sync
perl /data/home/bty565/empirical/inputs/popoolation2_1201/mpileup2sync.pl --input cage_1.mpileup --min-qual 20 --output cage_1.sync # CHANGE AS REQUIRED