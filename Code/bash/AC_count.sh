#Code to generate pileup files from BAM, match against all intronic regions, calculate AFs, and calculate ACs for input into the NB software

#Module load for cluster
module load samtools
module load julia

#Unzip pileup files
cd /data3/scratch/bty565/bams/cage_1/pileup
#echo 'Uncompressing pileup files...'
#for x in *.pileup.gz ; do gunzip "$x" ; done
#echo 'Extraction complete!'

#Grab all variants within intron regions
echo 'Grabbing intronic variants...'
mkdir allele_count
for x in *.pileup ; do 
  while read -r chr start end ; do
    awk -v chr="$chr" -v start="$start" -v end="$end" '
      $1 == chr && $2 >= start && $2 <= end {
        print
      }
    ' "$x"
  done < intron_pos.txt ;
done > allele_count/"$x"
echo 'intronic variants extracted!'

#Run ngsPool to calculate minor allele frequencies
echo 'Running ngsPool......'
cd allele_count
mkdir ngs_output
for x in *.pileup ; 
  do julia /data/home/bty565/empirical/ngsJulia/ngsPool/ngsPool.jl --fin "$x" --fout ngs_output/"$x".out.gz --nChroms 96 ;
  done
echo 'Allele frequencies calculated!'

#-----TODO-----
#This is the step in which the filtering step should be coded in! Make sure that is done by the notes in the notes file and then you can continue with the next steps for the code
#Check message related to file name length

#Extract column 11 (minor allele frequency) and output to new file
cd ngs_output/
mkdir ac_files
echo 'Extracting 11 for position and minor allele frequency......'
for x in *.out.gz ;
  do zcat "$x" | awk '{print$11}' > ac_files/"$x" ;
  done
cd ac_files/
#Remove saf_MLE string from the first line of each file. Removes the whole line 
for x in *.gz ;
  do sed -i '/saf_MLE/d' "$x" ;
  done
echo 'Output file created!'

#Read into each file and multiple AF by 96
for y in *.gz ;
  do while read -r x ; 
    do echo | awk -v af="$x" 'BEGIN {print af * 96}' ;
      done < "$y" > "$y"_minor.txt ;
  done
#TODO...test line of code above

#Read into each output from step above and substract from 96 and output to a new file
for y in *_minor.txt ;
  do while read -r x ;
    do echo | awk -v af="$x" 'BEGIN {print 96 - af}' ;
  done < "$y" > "$y"_major.txt ;
done

#Paste major and minor files together and convert tab into space
#TODO....

#Merges all text files together
#for f in *.csv ; do (cat "${f}"; echo) >> finalfile.txt; done
