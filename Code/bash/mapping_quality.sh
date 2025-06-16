# Script to read into Binary Alignment Files and extract mapping quality (MAPQ) and genomic position
# Run https://github.com/BaronKoy/MLBS/blob/main/Code/r/mapping.r script to plot results (csv input file will need to be changed accordingly

# Convert bams to sams
mkdir sams
for x in *.bam ; do samtools view "$x" > sams/*.sam ; done # Change location if running from cluster

# Extract chromosome 3R and output MAPQ and genomic position to new file
for x in *.sam ; do awk '$3=="3R"' "$x" | awk ‘{print$4,$5}’ > "$x"_region.txt ; done

# Extract 1000bp focal region
# TODO...debug so even if region string is missing from the file, the closest region is still found
mkdir output
for x in *_region.txt  ; do awk '/18520500/,/18521500/' "$x" > output/"$x"
cd output

# Add generation number to individual files
for x in S??0817_* ; do sed -i 's/$/,2/g' "$x" ; done
for x in S??0917_* ; do sed -i 's/$/,4/g' "$x" ; done
for x in S??1117_* ; do sed -i 's/$/,8/g' "$x" ; done
for x in S??0118_* ; do sed -i 's/$/,12/g' "$x" ; done
for x in S??0518_* ; do sed -i 's/$/,20/g' "$x" ; done
for x in S??0918_* ; do sed -i 's/$/,28/g' "$x" ; done
for x in S??0119_* ; do sed -i 's/$/,36/g' "$x" ; done
for x in S??0519_* ; do sed -i 's/$/,44/g' "$x" ; done
for x in S??1119_* ; do sed -i 's/$/,56/g' "$x" ; done

# Convert space between position and generation into comma
for x in *.txt ; do sed -i ‘s/ /,/g’ "$x" ; done