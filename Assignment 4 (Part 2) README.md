# My-repo-SLE777
ASSIGNMENT 4 (Part 2)- Biological Sequence Diversity Analysis


# Overview
In the second part of the project, coding DNA sequences (CDS) from Escherichia coli and Saprospirales bacteria are used to study biological sequence diversity.  As part of the jobs, you will need to download CDS data, count genes and total coding DNA, and look at sequence lengths, nucleotide and amino acid makeup, codon usage bias, and protein k-mer enrichment.  The goal is to use bioinformatics tools to look into how the genomes of two different types of bacteria are different.  The process is fully automated, well-marked, and can be used again and again in any R setting.

# Repository Structures
The project files are the same as in Part 1 so that everything is the same.  The raw sequence data are kept in data/, and all the outputs, like lists and plots, are kept in outputs/part2.  There are six main R scripts in the scripts/ folder. Each one is for one or more assignments.  The main report file report. Rmd takes these results and puts them all together into a single, organized file.  All scripts are made up of separate parts that can be run together or separately.

# Data Sources
suppressPackageStartupMessages({
  library("seqinr") # is a package designed to process and analyse sequence data.
  library("R.utils") # general utilities like zip and unzip

  library("R.utils")

URL="http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
download.file(URL,destfile="ecoli_cds.fa.gz")

list.files()

library("R.utils")
URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_58_collection/saprospirales_bacterium_gca_003448025/cds/Saprospirales_bacterium_gca_003448025.ASM344802v1.cds.all.fa.gz"
download.file(URL,destfile="sapro_cds.fa.gz")

list.files()

library("seqinr")
cds <- seqinr::read.fasta("ecoli_cds.fa")
str(head(cds))

cds <- seqinr::read.fasta("sapro_cds.fa")
str(head(cds))

# How to run
# 1. Counting coding sequences-
nA <- length(cdsA)
nB <- length(cdsB)

q1_table <- data.frame(
  Organism  = c(orgA_name, orgB_name),
  CDS_count = c(nA, nB),
  row.names = NULL
)
print(q1_table)

Optional: save
write.csv(q1_table, "Q1_cds_counts.csv", row.names = FALSE)


