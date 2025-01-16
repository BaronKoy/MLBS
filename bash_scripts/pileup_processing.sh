mkdir outputs
for x in S01* ; do zcat "$x" | sed -n -e '/^3R/p' > outputs/"$x" ; done
mv outputs
for x in S01* ; do awk '/18364273/,/18595586/' "$x"_com. ; done

