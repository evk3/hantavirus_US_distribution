#!/bin/bash -l
## Grid Engine Example Job Script  

# -- Begin SGE embedded arguments --
#$ -V
#Pass all environment variables to job
#$ -cwd
#Use current working directory

#$ -N DeNovo
# Name of script

#$ -j y
#Combine standard error and output files.

#$-q highmem.q
#Use the highmem.q queue, and not any other queue.

#$-pe smp 2
#Ask for a parallel environment for multi-threading.

# -- End SGE embedded arguments --

module load bwa/0.7.17
module load samtools/1.10
module load cutadapt/2.3
module load prinseq/0.20.3
module load picard/2.23.0
module load BEDTools/2.27.1
module load SPAdes/3.15.3

module load Python3/3.7
module load java/latest

# create temp directory for work on /scratch
mkdir -p /scicomp/scratch/evk3/US_hanta/
scratch='/scicomp/scratch/evk3/US_hanta/'

echo "Node number: " $HOSTNAME

SEEDFILE=./file_names_de_novo2.txt
file_num=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

OUTPUT_DIR=/scicomp/home/evk3/Diagnostics/Hantavirus/Distribution_paper_Nov2021/de_novo_assembly
echo "Output directory is: " $OUTPUT_DIR

echo "Patient number" $file_num
echo "SGE Value " $SGE_TASK_ID

sample_num=$(sed -r 's/\.\/raw_data\/(.*)_.*$/\1/g' <<< "$file_num")
echo "Sample number is: " $sample_num

L1_READ1=$file_num
L1_READ2=$(awk "NR==($SGE_TASK_ID + 1)" $SEEDFILE)

echo $L1_READ1
echo $L1_READ2

cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna /"$scratch"/
cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna.amb /"$scratch"/
cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna.ann /"$scratch"/
cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna.bwt /"$scratch"/
cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna.pac /"$scratch"/
cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna.sa /"$scratch"/
cp ~/Diagnostics/Human_GRCh38/GCA_000001405.15_GRCh38_top-level.fna.fai /"$scratch"/

#Unzip files for prinseq-lite
echo "Gunzipping now!"
gunzip -c $L1_READ1 > "$scratch"/"$sample_num"_R1_cutadapt.fastq
gunzip -c $L1_READ2 > "$scratch"/"$sample_num"_R2_cutadapt.fastq

#Remove low quality reads:
echo "starting printseq-lite"
prinseq-lite -fastq "$scratch"/"$sample_num"_R1_cutadapt.fastq -fastq2 "$scratch"/"$sample_num"_R2_cutadapt.fastq -min_qual_mean 25 -trim_qual_right 20 -min_len 50 -out_good "$scratch"/"$sample_num"_trimmed

#Map to human reference, use unmapped reads for SPAdes.
echo "Mapping to human genome."
bwa mem -t $NSLOTS "$scratch"/GCA_000001405.15_GRCh38_top-level.fna "$scratch"/"$sample_num"_trimmed_1.fastq "$scratch"/"$sample_num"_trimmed_2.fastq > "$scratch"/"$sample_num"_L.sam

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_L.bam "$scratch"/"$sample_num"_L.sam

echo "Removing reads mapped to human genome."
samtools view -b -f 4 "$scratch"/"$sample_num"_L.bam > "$scratch"/"$sample_num"_L-unmapped.bam

echo "Convert unmapped bam file to fastq files"
picard SamToFastq INPUT="$scratch"/"$sample_num"_L-unmapped.bam FASTQ="$scratch"/"$sample_num"_unaligned_1.fastq SECOND_END_FASTQ="$scratch"/"$sample_num"_unaligned_2.fastq

#De novo assembling:
echo "de novo assembling!"
spades.py -k auto -t $NSLOTS -o "$scratch"/"$sample_num"_contigs -1 "$scratch"/"$sample_num"_unaligned_1.fastq -2 "$scratch"/"$sample_num"_unaligned_2.fastq

#Quast for metrics:
echo "Running quast for de novo assembly metrics:"
module unload Python3/3.7
module load Python3/3.4
quast.py --output-dir "$scratch"/"$sample_num"_contigs/quast "$scratch"/"$sample_num"_contigs/contigs.fasta
module unload Python3/3.4

cp "$scratch"/"$sample_num"_contigs/contigs.fasta "$OUTPUT_DIR"/"$sample_num"_contigs.fasta

module load Python3/3.6.1
#Identify contigs by nucleotide blast:
echo "Identifying contigs by nucleotide blast"
python3 blast_python.py -i "$scratch"/"$sample_num"_contigs/contigs.fasta -o "$scratch"/"$sample_num"
#python3 parse xml and return patient_num_blast_results.txt


cp  "$scratch"/"$sample_num".xml "$OUTPUT_DIR"/

cp  "$scratch"/"$sample_num"_blast_hits.csv "$OUTPUT_DIR"/

cp  "$scratch"/"$sample_num"_contigs/quast/report.txt "$OUTPUT_DIR"/"$sample_num"_report.txt


#clean up scratch directory
#rm -rf /scratch/evk3

module unload bwa/0.7.17
module unload samtools/1.10
module unload cutadapt/2.3
module unload prinseq/0.20.3
module unload picard/2.23.0
module unload BEDTools/2.27.1
module unload SPAdes/3.14.0
module unload Python3/3.7
module unload java/latest

echo "Script finish"
