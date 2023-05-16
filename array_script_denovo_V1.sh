#!/bin/bash -l

#The purpose of this wrapper script is to submit
#the deNovo_and_identify_V1.sh script using all files in file_names_de_novo.txt


files=$(awk 'END{print NR}' ./file_names_de_novo.txt)
echo $files

#Submit script from 1 to the total number of R1/R2.fastq.gz files, skip every other line when submitting (to account for the R1/R2 read pairing).
qsub -t 1-$files:2 \
	-N Identify \
	./deNovo_and_identify_V1.sh \

echo "Wrapper script done"
