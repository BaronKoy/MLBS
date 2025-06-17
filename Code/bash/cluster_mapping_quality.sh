# Script to read into Binary Alignment Files and extract mapping quality (MAPQ) and genomic position
# Run https://github.com/BaronKoy/MLBS/blob/main/Code/r/mapping.r script to plot results (csv input file will need to be changed accordingly

# Modules
module load samtools

# Convert bams to sams
wd=/data3/scratch/bty565/bams/cage_2
mkdir $wd/sams
wds=$wd/sams
for x in $wd/*.bam ; do samtools view "$x" > $wd.sam ; done

# Extract chromosome 3R and output MAPQ and genomic position to new file
for x in $wd/*.sam ; do awk '$3=="3R"' "$x" | awk ‘{print$4,$5}’ > "$x"_region.txt ; done

# Extract 1000bp focal region
# TODO...debug so even if region string is missing from the file, the closest region is still found
mkdir $wds/output
wdso=$wds/output
for x in $wd/*_region.txt  ; do awk '/18520500/,/18521500/' "$x" > "$x"_final.txt ; done
for x in $wd/*_final.txt ; do mv "$x" $wd/sams/output/"$x" ; done

# Add generation number to individual files
cd $wdso
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
for x in *.final.txt ; do sed -i ‘s/ /,/g’ "$x" ; done
