### ReproHackathon Project

# Group Members: Anakim GUALDONI - Georges HANNA - Noémie BOZIER - Camille DE AMORIM

# 1- Objective :
A Snakemake workflow that allows users to reproduce the results of the scientific article ( https://www.nature.com/articles/s41467-020-15966-7).

# 2- Usage :
.................

# 3- How does it work ?

First, the pipeline downloads all FASTQ files studied in the article using fasterq-dump (version 3.1.1).

Next, the FASTQ files are trimmed using Trim Galore (version ___). During this step, bases with a quality score below 20 will be trimmed from the ends of the reads, and only reads with a minimum length of 25 bases will be retained.

The pipeline also downloads the reference genome in FASTA format and the genome annotation in GFF format. 
Finally, it ensures the indexing of the reference genome using Bowtie (version ___).