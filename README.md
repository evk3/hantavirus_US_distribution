# hantavirus_US_distribution

The site contains the code used in this manuscript: "Hantavirus Prevalence and Viral Strain Distribution in the US, 2008-2020"

## PDF Parser
This script parsed the results of CDC pdf reports sent to submitters (state public health labs and requesting doctors) to compile hantavirus results from 2008-2020.

## Download all available hantavirus genomes.
This script is used to download all available hantavirus genomes.  The script can be easily modified to download other viral genomes, too!

## Bioinformatics
This code was used to build unknown and diverse hantavirus genomes.  Reads were quality trimmed, mapped to a human reference and unmapped reads were de novo assembled with SPaDes.  Contigs were blasted and the resulting hantavirus hits were used as reference sequences for read mapping.  Initial scaffold genomes were built using both contigs and reads, and genomes were iteratively assembled using scaffold genomes and reads.  Final genomes were parsed by eye to identify any frameshifts or indels.  These troublemaker genomes were re-assembled using the Geneious mapper and reads (but not using contigs, since these introduced most of the errors).
