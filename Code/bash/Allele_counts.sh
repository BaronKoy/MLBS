#Code to generate allele frequencies for all non-coding variants and calculate the AC required for input into the NB software
#This assumes the input files are pileups which contain all non-coding variants

#TODO...currently the script needs to be run from the same working directory the pileup files reside in

#Unzip pileup files
echo 'Uncompressing pileup files...'
for x in *.mpileup.gz ; do gunzip "$x" ; done
echo 'Extraction complete!'

#Match raw pileup files against non-coding variant file, output matches to separate file
#Run in the directory containing the raw pileup files
echo 'Matching non-coding variants against uncompressed pileups'
for x in *.mpileup ; do grep -wFf ../non_coding.txt "$x" > ../cage_1/"$x" ; done

#Change and create working directories
echo 'Changing into working directories'
cd ../cage_1
mkdir ngs_output

#Run ngsPool to calculate minor allele frequencies
echo 'Running ngsPool......'
for x in *.pileup ; do julia ~/software/ngsJulia/ngsPool/ngsPool.jl --fin "$x" --fout ngs_output/"$x".out.gz --nChroms 96 ; done
echo Allele frequencies calculated!

#-----TODO-----
#This is the step in which the filtering step should be coded in! Make sure that is done by the notes in the notes file and then you can continue with the next steps for the code

#Extract column 11 (minor allele frequency) and output to new file
cd ngs_output/
mkdir ac_files
echo 'Extracting 11 for position and minor allele frequency......'
for x in *.out.gz ; do zcat "$x" | awk '{print$11}' > ac_files/"$x" ; done
cd ac_files/
#Remove saf_MLE string from the first line of each file. Removes the whole line 
for x in *.gz ; do sed -i 's/saf_MLE/d' "$x" ; done
echo 'Output file created!'

#Read into each file and multiple AF by 96
for y in *.gz ; do while read -r x ; do echo | awk -v af="$x" 'BEGIN {print af * 96}' ; done < "$y" > "$y"_minor.txt ; done
#TODO...test line of code above

#Read into each output from step above and substract from 96 and output to a new file
for y in *_minor.txt ; do while read -r x ; do echo | awk -v af="$x" 'BEGIN {print 96 - af}' ; done < "$y" > "$y"_major.txt ; done

#Paste major and minor files together and convert tab into space
#TODO....

#Merges all text files together
for f in *.txt; do (cat "${f}"; echo) >> finalfile.txt; done

#Code blocks below used for all Cages
"""
# Extract colum 2 (Position) & column 11 (minor allele frequency) and output to new file
echo 'Extracting columns 2 & 11 for position and minor allele frequency......'
cd ngs_output/
for x in *.out.gz ; do zcat  | awk '{print$2,$11}' > plot_data/ ; done
echo 'Output file created!'

# Converting to csv
echo 'Converting to .csv......'
cd plot_data
for x in *.out.gz ; do sed -i 's/ /,/g'  ; done
echo 'Converted to csv file!'

# Add generation numbers to individual files
echo 'Adding generations to input files......'
for x in S??0817_* ; do sed -i 's/$/,2/g'  ; done
for x in S??0917_* ; do sed -i 's/$/,4/g'  ; done
for x in S??1117_* ; do sed -i 's/$/,8/g'  ; done
for x in S??0118_* ; do sed -i 's/$/,12/g'  ; done
for x in S??0518_* ; do sed -i 's/$/,20/g'  ; done
for x in S??0918_* ; do sed -i 's/$/,28/g'  ; done
for x in S??0119_* ; do sed -i 's/$/,36/g'  ; done
for x in S??0519_* ; do sed -i 's/$/,44/g'  ; done
for x in S??1119_* ; do sed -i 's/$/,56/g'  ; done
echo 'Generation numbers added to output files!'

# Removing non-variant lines
echo 'Removing non-variant lines......'
for x in *.out.gz ; do sed -i '/,0.0,/d'  ; done
for x in *.out.gz ; do sed -i '/position,saf_MLE/d'  ; done
echo 'Non-variant lines removed!'

# Combine into a single plot file
#TODO...concatenates final files into plot data file and add column headers for plot script..."""
