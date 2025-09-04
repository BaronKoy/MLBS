# Code to generate and cut pileup files to only include focal region 3R:18520501 - 18521500
# Retrives position, minor allele frequency, & adds generation to each file
# Produces final outfile file. Run Code/r/allele_frequency.r on final file to generate frequency plot

# Modules for cluster run
module load samtools
module load julia

# BAM to pileup
cd /data3/scratch/bty565/cages/cage_1/bams
mkdir pileup
for x in *.bam ; do samtools mpileup "$x" --max-depth 0 --min-BQ 0 --fasta-ref /data/home/bty565/empirical/inputs/GCA_000001215.4_full.fna -o pileup/"$x".pileup.gz ; done # change FASTA depending on location

# Extract lines for focal region and output to new file
cd pileup
mkdir frequency
echo 'Cutting mpileup files for focal region 3R:28520501-18521500......'
for x in *.gz ; do zcat "$x" | sed -rn '/^3R/p' | awk '/18520501/,/18521500/' > frequency/"$x".pileup ; done
echo Region cut complete!

# Run ngsPool to calculate minor allele frequencies
echo 'Running ngsPool......'
cd frequency
mkdir ngs_output
for x in *.pileup ; do julia /data/home/bty565/empirical/ngsJulia/ngsPool/ngsPool.jl --fin "$x" --fout ngs_output/"$x".out.gz --nChroms 96 ; done
echo Allele frequencies calculated!

# Extract colum 2 (Position) & column 11 (minor allele frequency) and output to new file
echo 'Extracting columns 2 & 11 for position and minor allele frequency......'
cd ngs_output/
mkdir plot_data
for x in *.out.gz ; do zcat "$x" | awk '{print$2,$11}' > plot_data/"$x" ; done
echo 'Output file created!'

# Converting to csv
echo 'Converting to .csv......'
cd plot_data
for x in *.out.gz ; do sed -i 's/ /,/g' "$x" ; done
echo 'Converted to csv file!'

# Add generation numbers to individual files
echo 'Adding generations to input files......'
for x in S??0817_* ; do sed -i 's/$/,2/g' "$x" ; done
for x in S??0917_* ; do sed -i 's/$/,4/g' "$x" ; done
for x in S??1117_* ; do sed -i 's/$/,8/g' "$x" ; done
for x in S??0118_* ; do sed -i 's/$/,12/g' "$x" ; done
for x in S??0518_* ; do sed -i 's/$/,20/g' "$x" ; done
for x in S??0918_* ; do sed -i 's/$/,28/g' "$x" ; done
for x in S??0119_* ; do sed -i 's/$/,36/g' "$x" ; done
for x in S??0519_* ; do sed -i 's/$/,44/g' "$x" ; done
for x in S??1119_* ; do sed -i 's/$/,56/g' "$x" ; done
echo 'Generation numbers added to output files!'

# Removing non-variant lines
echo 'Removing non-variant lines......'
for x in *.out.gz ; do sed -i '/,0.0,/d' "$x" ; done
for x in *.out.gz ; do sed -i '/position,saf_MLE/d' "$x" ; done
echo 'Non-variant lines removed!'

# Combine into a single plot file
echo 'Combining files into final dataset......'

# Define output filename
output_file="final_plot_data.csv"

# Add header row
echo "Position,Allele_frequency,Generation" > "$output_file"

# Concatenate files in correct generation order, sorting by position
for gen in 2 4 8 12 20 28 36 44 56 ; do
    grep -h ",$gen$" *.out.gz | sort -t',' -k1,1n >> "$output_file"
done

echo "Final dataset created: $output_file"
