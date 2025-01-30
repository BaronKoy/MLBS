# Cut original pileup files to only include focal region 3R:18520501 - 18521500
# TODO...currently the script needs to be run from the same working directory the data files reside in
# For the love of god fix the newline shell syntax in this

echo Current working directory: /home/baron/Documents/PhD/ML_balancingselection_data/abi_raw_data
echo Cutting mpileup files for focal region 3R:28520501-18521500...
for x in *.gz ; do zcat "$x" | sed -rn '/^3R/p' | awk '/18520501/,/18521500/' > region/"$x" ; done
echo Region cut complete!

# Run ngsPool for allele frequencies
echo Moving into region directory
cd region
echo Current working directory: /home/baron/Documents/PhD/ML_balancingselection_data/abi_raw_data/region
echo Running ngsPool...
for x in *.mpileup ; do julia ~/software/ngsJulia/ngsPool/ngsPool.jl --fin "$x" --fout ngspool_out/"$x" --nChroms 96
echo ngsPool process complete!