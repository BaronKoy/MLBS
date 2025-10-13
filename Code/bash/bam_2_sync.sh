#!/usr/bin/env bash

module load samtools
module load bedtools2

# change to working directory - specific cage - TODO when integrating into pipeline this should be an input variable
cd /data3/scratch/bty565/cages/cage_1 # Change as required

# Filter bams to only include neutral SNPs
for x in *.bam ; do samtools index "$x" ; done
for x in *.bam ; do bedtools intersect -abam "$x" -b /data/home/bty565/empirical/inputs/introns.bed > "$x"_intersect.bam ; done

# Collect files
files=(S01*.bam_intersect.bam) # Change prefix as required

declare -A date_map

for file in "${files[@]}"; do
  # Extract last 4 digits after the initial prefix
  mm_yy=$(echo "$file" | sed -E 's/^S01([0-9]{4}).*/\1/') # Change prefix as required
  
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
  mv "$file" "${count}.bam"
done