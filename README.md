# hantavirus_US_distribution

The site contains the code used in this manuscript: "Hantavirus Prevalence and Viral Strain Distribution in the US, 2008-2020"

## PDF Parser
The pdf_parsing_script.ipynb is an iPython Notebook that parses the CDC pdf reports sent to submitters (state public health labs and requesting doctors) to compile hantavirus results from 2008-2020.

## Download all available hantavirus genomes.
The DL_Hanta_genomes.ipynb script is used to download all available hantavirus genomes.  The script can be easily modified to download other viral genomes, too!  Run cells sequentially to download genomes from NCBI/Genbank and then parse these results into segment-specific Excel files.  The last cell contains some bonus, extra-handy code for converting a fasta file into a relaxed phylip file.  

## Pan-Hantavirus probe development and design.
All available hantavirus genomes were downloaded and duplicate sequences were removed by hand.  Hantavirus genomes were parsed using https://github.com/evk3/Nipah_phylogenetics to produce 5,325 probes that are 80bp long and  5' biotinylated.  Probe sequences can be found in "NGS submission Form_Whitmer_hanta.xlsx."  I order probes from Twist Biosciences <https://www.twistbioscience.com/>.

## Bioinformatics
Sharing the code that was used to build unknown and diverse hantavirus genomes.  

### Step 1: Identify the closest reference sequence.
Reads were quality trimmed, mapped to a human reference and unmapped reads were de novo assembled with SPaDes.  Contigs were blasted and the resulting hantavirus hits were used as reference sequences for read mapping.

-Submit script to SGE scheduler node: "bash array_script.sh"

-array script reads the file names in file_names.txt and submits each job, but skips every other line in file_names.txt to account for the R1/R1.fastq.gz files.

-file_name.txt format:  ./raw_data/sample1_S1_L001_R1_001.fastq.gz
                        ./raw_data/sample1_S1_L001_R2_001.fastq.gz
                        ./raw_data/sample2_S2_L001_R1_001.fastq.gz
                        ./raw_data/sample2_S2_L001_R2_001.fastq.gz
                        ....etc...
                        
 -deNovo_and_identify.sh runs for each job with input files: sample1_S1_L001_R1_001.fastq.gz and sample1_S1_L001_R2_001.fastq.gz.
        -data are trimmed and de novo assembled.
        
 -blast_python.py performs the blasting of contigs and returns the results in an *.xml and *.txt format.

### Step 2: Build genomes using the closest matching reference sequence.
Initial scaffold genomes were built using both contigs and reads, and genomes were iteratively assembled using scaffold genomes and reads.  Final genomes were parsed by eye to identify any frameshifts or indels.  These troublemaker genomes were re-assembled using the Geneious mapper and reads (but not using contigs, since these introduced most of the errors).
