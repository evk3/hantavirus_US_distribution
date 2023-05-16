#!/bin/bash -l
## Grid Engine Example Job Script  

# -- Begin SGE embedded arguments --
#$ -V
#Pass all environment variables to job
#$ -cwd
#Use current working directory

#$ -j y
#Combine standard error and output files.

#$-q short.q
#Use the short.q queue, and not any other queue.

#$-pe smp 2
#Ask for a parallel environment for multi-threading.

#$-l h=!'node228.hpc.biotech.cdc.gov'&!'node229.hpc.biotech.cdc.gov'&!'node230.hpc.biotech.cdc.gov'&!'node46.hpc.biotech.cdc.gov'

# -- End SGE embedded arguments --

module load bwa/0.7.17
module load samtools/1.14
module load picard/2.23.0
module load BEDTools/2.27.1
module load cutadapt/2.3
module load prinseq/0.20.3
module load htslib/1.9
module load bowtie2/2.3.5.1
module load gatk/4.1.7.0
module load gcc/9.2.0
module load htslib/1.10

# create temp directory for work on /scratch

mkdir -p /scicomp/scratch/evk3/US_hanta/
scratch='/scicomp/scratch/evk3/US_hanta'

echo "Hostname:" $HOSTNAME
echo "SGE Value " $SGE_TASK_ID

OUTPUT_PATH=/scicomp/home/evk3/Diagnostics/Hantavirus/Distribution_paper_Nov2021/mapping_output/
echo "Output path: " $OUTPUT_PATH

REFERENCE_PATH=$(awk "NR==($SGE_TASK_ID)" ./references_read-contig_mapping_hanta.txt)
cp $REFERENCE_PATH $scratch/
REFERENCE=$(sed -r 's/\.\/references\/(.*$)/\1/g' <<< "$REFERENCE_PATH")

echo "Reference: " $REFERENCE

SEEDFILE=./file_names_hanta.txt
file_num=$(awk "NR==$SGE_TASK_ID" $SEEDFILE)

echo "First File" $file_num

sample_num=$(sed -r 's/\.\/raw_data\/(.*)_.*_.*_.*_.*$/\1/g' <<< "$file_num")
echo "Sample number is: " $sample_num

L1_READ1=$file_num
L1_READ2=$(awk "NR==($SGE_TASK_ID + 1)" $SEEDFILE)

echo $L1_READ1
echo $L1_READ2

cp ./de_novo_assembly/"$sample_num"_contigs.fasta /"$scratch"/
CONTIGS=/"$scratch"/"$sample_num"_contigs.fasta


#Remove Illumina TruSeq adaptors from PE reads:
echo "Starting cutadapt"
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
         -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         --cores=$NSLOTS \
         -m 1 \
         -o "$scratch"/"$sample_num"_R1_cutadapt.fastq.gz -p "$scratch"/"$sample_num"_R2_cutadapt.fastq.gz \
         $L1_READ1 $L1_READ2 

echo "Gunzipping now!"
gunzip -c "$scratch"/"$sample_num"_R1_cutadapt.fastq.gz > "$scratch"/"$sample_num"_R1_cutadapt.fastq
gunzip -c "$scratch"/"$sample_num"_R2_cutadapt.fastq.gz > "$scratch"/"$sample_num"_R2_cutadapt.fastq

#Remove low quality reads:
echo "starting printseq-lite"
prinseq-lite -fastq "$scratch"/"$sample_num"_R1_cutadapt.fastq -fastq2 "$scratch"/"$sample_num"_R2_cutadapt.fastq -min_qual_mean 25 -trim_qual_right 20 -min_len 50 -out_good "$scratch"/"$sample_num"_trimmed

echo "Indexing reference sequence using bowtie2"
bwa index "$scratch"/"$REFERENCE" -p "$scratch"/"$REFERENCE"

echo "Mapping reads to reference genome"
bwa mem -t $NSLOTS "$scratch"/"$REFERENCE" "$scratch"/"$sample_num"_trimmed_1.fastq "$scratch"/"$sample_num"_trimmed_2.fastq > "$scratch"/"$sample_num"_reads.sam

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_reads.bam "$scratch"/"$sample_num"_reads.sam

#******************************************************************************************************
#If contigs file exists, map contigs.
#else, just use reads-based mapping.
if test -f "$CONTIGS"
then
	echo "Mapping contigs to reference genome"
	bwa mem -t $NSLOTS "$scratch"/"$REFERENCE" "$scratch"/"$sample_num"_contigs.fasta  > "$scratch"/"$sample_num"_contigs.sam
	
	echo "Starting samtools - convert SAM to BAM"
	samtools view -S -b -o "$scratch"/"$sample_num"_contigs.bam "$scratch"/"$sample_num"_contigs.sam

	echo "Merging read and contig bam files"
	samtools merge -f "$scratch"/"$sample_num"_merge.bam "$scratch"/"$sample_num"_reads.bam "$scratch"/"$sample_num"_contigs.bam
else
	echo "In else statement, no contings to map, mapping reads only!"
	mv "$scratch"/"$sample_num"_reads.bam "$scratch"/"$sample_num"_merge.bam
fi

# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing Ebo genome with samtools"
samtools faidx "$scratch"/"$REFERENCE"


echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_merge.bam -o "$scratch"/"$sample_num"_merge.sorted.bam

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_merge.sorted.bam


JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp

#Copy only mapped bases, and save:
samtools view -b -F 4 "$scratch"/"$sample_num"_merge.sorted.bam > "$scratch"/"$sample_num"_merge-mapped.sorted.bam


#Make intermediate1 fasta:
echo "Making intermediate1 fasta!"
bcftools mpileup -A -d 6000000 -B -Q 0 -Ov -f "$scratch"/"$REFERENCE" "$scratch"/"$sample_num"_merge-mapped.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o "$scratch"/"$sample_num"_merge-mapped.sorted.bam.vcf

gatk IndexFeatureFile --input "$scratch"/"$sample_num"_merge-mapped.sorted.bam.vcf

gatk CreateSequenceDictionary -R "$scratch"/"$REFERENCE"

gatk FastaAlternateReferenceMaker -R "$scratch"/"$REFERENCE" -O "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta -V "$scratch"/"$sample_num"_merge-mapped.sorted.bam.vcf

#**********************************************************************************************************

echo "Performing mapping to Self #1"

echo "Indexing reference sequence using bowtie2"
bowtie2-build-s "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta "$scratch"/"$sample_num"_intermediate1_myIndex 

echo "Mapping reads to reference genome"
bowtie2-align-s -I 0 -X 800 -p 32 --sensitive -q -x "$scratch"/"$sample_num"_intermediate1_myIndex -1 "$scratch"/"$sample_num"_trimmed_1.fastq -2 "$scratch"/"$sample_num"_trimmed_2.fastq -S "$scratch"/"$sample_num"_intermediate1.sam

# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing inetermediate1 genome with samtools"
samtools faidx "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_intermediate1.bam "$scratch"/"$sample_num"_intermediate1.sam

echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_intermediate1.bam -o "$scratch"/"$sample_num"_intermediate1.sorted.bam

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_intermediate1.sorted.bam


JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp

#Copy only mapped bases, and save:
samtools view -b -F 4 "$scratch"/"$sample_num"_intermediate1.sorted.bam > "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam


#Make intermediate1 fasta:
echo "Making intermediate1 fasta!"
bcftools mpileup -A -d 6000000 -B -Q 0 -Ov -f "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam | bcftools call --ploidy 1 -mv -Ov -o "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam.vcf

gatk IndexFeatureFile --input "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam.vcf

gatk CreateSequenceDictionary -R "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta

gatk FastaAlternateReferenceMaker -R "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta -O "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta -V "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam.vcf

#**********************************************************************************************************

echo "Performing mapping to Self #2"

echo "Indexing reference sequence using bowtie2"
bowtie2-build-s "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$scratch"/"$sample_num"_intermediate2_myIndex 

echo "Mapping reads to reference genome"
bowtie2-align-s -I 0 -X 800 -p 32 --sensitive -q -x "$scratch"/"$sample_num"_intermediate2_myIndex -1 "$scratch"/"$sample_num"_trimmed_1.fastq -2 "$scratch"/"$sample_num"_trimmed_2.fastq -S "$scratch"/"$sample_num"_intermediate2.sam

# Previously generated samtools index of reference genome.  Generates *.fai file and only need to do 1X.
# This step may not be necessary with the samtool pipeline below:
echo "Indexing inetermediate1 genome with samtools"
samtools faidx "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta

echo "Starting samtools - convert SAM to BAM"
samtools view -S -b -o "$scratch"/"$sample_num"_intermediate2.bam "$scratch"/"$sample_num"_intermediate2.sam

echo "Starting samtools sort BAM file"
samtools sort -@ $NSLOTS "$scratch"/"$sample_num"_intermediate2.bam -o "$scratch"/"$sample_num"_intermediate2.sorted.bam

echo "Starting samtools index BAM file"
samtools index "$scratch"/"$sample_num"_intermediate2.sorted.bam


JAVA_OPTS='-Xmx50g'
TMP_DIR=/tmp

#Copy only mapped bases, and save:
samtools view -b -F 4 "$scratch"/"$sample_num"_intermediate2.sorted.bam > "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam

samtools index "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam

#Make intermediate1 fasta:
echo "Making final fasta!"

echo "Making consensus fasta!"

samtools mpileup -r 1 -A -aa -d 6000000 -B -Q 0 -f "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam | /scicomp/home-pure/evk3/setup/ivar_1.3.1/src/ivar consensus -p "$scratch"/"$sample_num".L.consensus -m 10 -n N

samtools mpileup -r 2 -A -aa -d 6000000 -B -Q 0 -f "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam | /scicomp/home-pure/evk3/setup/ivar_1.3.1/src/ivar consensus -p "$scratch"/"$sample_num".M.consensus -m 10 -n N

samtools mpileup -r 3 -A -aa -d 6000000 -B -Q 0 -f "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam | /scicomp/home-pure/evk3/setup/ivar_1.3.1/src/ivar consensus -p "$scratch"/"$sample_num".S.consensus -m 10 -n N


#**********************************************************************************************************


#copy results from node /scratch/evk3/ebo back to home dir
cp "$scratch"/"$sample_num"_intermediate1_no_contigs.fasta "$OUTPUT_PATH"
cp "$scratch"/"$sample_num"_intermediate2_no_contigs.fasta "$OUTPUT_PATH"
cp "$scratch"/"$sample_num"_merge-mapped.sorted.bam "$OUTPUT_PATH"
cp "$scratch"/"$sample_num"_intermediate1-mapped.sorted.bam "$OUTPUT_PATH"
cp "$scratch"/"$sample_num"_intermediate2-mapped.sorted.bam "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".L.consensus.fa "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".L.consensus.qual.txt "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".M.consensus.fa "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".M.consensus.qual.txt "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".S.consensus.fa "$OUTPUT_PATH"
cp "$scratch"/"$sample_num".S.consensus.qual.txt "$OUTPUT_PATH"

#clean up scratch directory
#rm -rf /scratch/evk3



module unload bwa/0.7.17
module unload samtools/1.9
module unload picard/2.21.1
module unload BEDTools/2.27.1
module unload cutadapt/2.3
module unload prinseq/0.20.3
module unload htslib/1.9
module unload bowtie2/2.3.5.1
module unload gatk/4.1.7.0
module unload htslib/1.10
module unload gcc/9.2.0

echo "Script finish"

