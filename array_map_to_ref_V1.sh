#!/bin/bash -l

#The purpose of this wrapper script is to submit
#the read_contig_iterative_map_V1.sh script using all files in file_names_hanta.txt


files=$(awk 'END{print NR}' ./file_names_hanta.txt)
echo $files

#Submit script from file 1 to total length, but skip every other file name (to account for files in the R1/R2 format).
qsub -t 1-$files:2 \
	-N Hanta \
	./read_contig_iterative_map_V1.sh \

echo "Wrapper script done"
