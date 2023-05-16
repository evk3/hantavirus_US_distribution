#!/bin/usr
#This is written in python3.

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
import getopt, sys
import os
from Bio import SeqIO
import re
import pandas as pd
from os import listdir

#get cwd:
current_directory=os.getcwd()

#read commandLine arguments, first
fullCmdArguments = sys.argv

# - further arguments
argumentList = fullCmdArguments[1:]

#print(argumentList)
#Usage: python3 blast_sigma_fasta_V2.py -i $file_num -o $file_string
#output: ['-i', 'sigma_16_mut0.fasta', '-o', 'sigma_16_mut0']

blastn_cline = NcbiblastnCommandline(cmd="/apps/x86_64/ncbi/ncbi-blast-2.9.0+/bin/blastn", query=argumentList[1], db="/scicomp/reference/ncbi-blast/old/v5/nt_v5/nt_v5", outfmt=5, out="%s.xml" % (argumentList[3]), task="blastn", max_target_seqs="1")
blastn_cline
print(blastn_cline)
stdout, stderr = blastn_cline()

print(stdout)
print(stderr)
print("Done blasting!")

print("Parsing and writing blast results")

blast_results = SearchIO.parse("%s.xml" % (argumentList[3]), "blast-xml")

df_columns=['contig_name', 'contig_length', 'accession_num', 'blast_hit']
df_Compiled_Data = pd.DataFrame(columns=df_columns)

total_blast_results=0
for blast_record in blast_results:
    seq_length=blast_record.seq_len
    sequence_name=blast_record.id
    
    for hit in blast_record:
        description=hit.description
        accession=hit.accession
        
        temp_df = pd.DataFrame({'contig_name':[sequence_name],
                                'contig_length':[seq_length], 
                                'accession_num':[accession],
                                'blast_hit':[description]})

        df_Compiled_Data=df_Compiled_Data.append(temp_df)
        total_blast_results = total_blast_results + 1
    
    if len(blast_record) == 0:
        #print("Found contig with no blast hit!")
        temp_df = pd.DataFrame({'contig_name':[sequence_name],
                                'contig_length':[seq_length], 
                                'accession_num':["XXXXXX"],
                                'blast_hit':["No blast hit found"]})

        df_Compiled_Data=df_Compiled_Data.append(temp_df)
        
#Print number of hits identified / total contigs made.
print("Total contigs generated: %s" % len(df_Compiled_Data))
print("Total number of identified contigs: %s" % total_blast_results)
print("%s contigs not identified" % (len(df_Compiled_Data) - total_blast_results))

#Write contents of data frame to text file.        
df_Compiled_Data.to_csv("%s_blast_hits.csv" % (argumentList[3]), sep='\t', index=True, header=True)

print("Wrote blast results to csv")

