# Cut original pileup files to only include focal region 3R:18520501 - 18521500
# TODO...currently the script needs to be run from the same working directory the data files reside in
# For the love of god fix the newline shell syntax in this

echo 'Cutting mpileup files for focal region 3R:28520501-18521500......'
for x in *.gz ; do zcat "$x" | sed -rn '/^3R/p' | awk '/18520501/,/18521500/' > "$x".pileup ; done
echo Region cut complete!

# Run ngsPool for allele frequencies
echo 'Running ngsPool......'
for x in *.pileup ; do julia ~/software/ngsJulia/ngsPool/ngsPool.jl --fin "$x" --fout ngs_output/"$x".out.gz --nChroms 96 ; done
echo ngsPool process complete!

# Extract colum 2 (Position) & column 11 (minor allele frequency)
echo 'Extracting columns 2 & 11 for position and minor allele frequency......'
cd ngs_output/
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

# Removing non-variant lines
echo 'Removing non-variant lines......'
for x in *.out.gz ; do sed -i '/,0.0,/d' "$x" ; done
for x in *.out.gz ; do sed -i '/position,saf_MLE/d' "$x" ; done
echo 'Non-variant lines removed!'

# Combine into a single plot file
#TODO...