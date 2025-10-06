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

# 2. Total coding DNA length-
Robust total-coding-nt calculator that works whether cdsA/cdsB are DNAStringSet, character, or list
get_total_nt <- function(x) {
  if (inherits(x, "XStringSet")) return(sum(width(x)))
  if (inherits(x, "XString"))    return(as.integer(width(x)))
  if (is.character(x))           return(sum(nchar(x)))
  if (is.list(x))                return(sum(nchar(unlist(lapply(x, as.character), use.names = FALSE))))
  stop("Unsupported type for cds object")
}

 Compute totals and make the table (expects orgA_name, orgB_name to be defined)
total_nt_A <- get_total_nt(cdsA)
total_nt_B <- get_total_nt(cdsB)

q2_table <- data.frame(
  Organism        = c(orgA_name, orgB_name),
  Total_coding_nt = c(total_nt_A, total_nt_B),
  stringsAsFactors = FALSE
)

print(q2_table)

# 3. CDS length distribution-
Calculate CDS lengths for each gene
cds_len_A <- vapply(cdsA, length, integer(1))
cds_len_B <- vapply(cdsB, length, integer(1))

 Summary stats
q3_stats <- data.frame(
  Organism  = c(orgA_name, orgB_name),
  Mean_nt   = c(mean(cds_len_A), mean(cds_len_B)),
  Median_nt = c(median(cds_len_A), median(cds_len_B)),
  SD_nt     = c(sd(cds_len_A), sd(cds_len_B)),
  Min_nt    = c(min(cds_len_A), min(cds_len_B)),
  Max_nt    = c(max(cds_len_A), max(cds_len_B))
)

print(q3_stats)

 Boxplot (on screen)
boxplot(list(E_coli = cds_len_A, Saprospirales = cds_len_B),
        ylab = "CDS length (nt)", 
        main = "CDS length distribution", 
        outline = FALSE)
grid()

 Optional: save figure
 png("Q3_cds_length_boxplot.png", 1100, 700)
 boxplot(list(E_coli = cds_len_A, Saprospirales = cds_len_B),
         ylab="CDS length (nt)", 
         main="CDS length distribution", 
         outline=FALSE)
 grid()
 dev.off()

 # 4. Base and amino-acid frequencies-
 helpers you need to add once (put above your code) ---
 Collapse a list of CDS (each is a vector of bases) into one long vector
collapse_dna <- function(cds_list) {
  unlist(cds_list, use.names = FALSE)
}

 Translate a CDS to amino acids safely
 (returns a vector of single-letter amino acids; NA if translation fails)
safe_translate <- function(dna_vec) {
  out <- tryCatch(seqinr::translate(dna_vec, numcode = 1),
                  error = function(e) NA_character_)
  out
}


 Nucleotide frequencies (pooled across all CDS)
valid_bases <- c("A","C","G","T")
count_nt <- function(vec) {
  tbl <- table(factor(toupper(vec), levels = valid_bases))
  as.numeric(tbl)
}

dnaA_all <- collapse_dna(cdsA)
dnaB_all <- collapse_dna(cdsB)

ntA <- count_nt(dnaA_all);  ntA_freq <- ntA / sum(ntA)
ntB <- count_nt(dnaB_all);  ntB_freq <- ntB / sum(ntB)

q4_nt <- data.frame(
  Base = valid_bases,
  Freq_E_coli = ntA_freq,
  Freq_Saprospirales = ntB_freq
)
print(q4_nt)

 Barplot (nucleotides) with colour
barplot(rbind(q4_nt$Freq_E_coli, q4_nt$Freq_Saprospirales),
        beside = TRUE,
        names.arg = q4_nt$Base,
        ylab = "Frequency",
        main = "Nucleotide frequencies (CDS total)",
        col = c("skyblue", "lightpink"))
legend("topright", legend = c("E. coli","Saprospirales"), bty = "n",
       fill = c("skyblue", "lightpink"))
grid()

 Amino-acid frequencies (translate each CDS, concatenate)
protA <- unlist(lapply(cdsA, safe_translate), use.names = FALSE)
protB <- unlist(lapply(cdsB, safe_translate), use.names = FALSE)
protA <- protA[!is.na(protA)]
protB <- protB[!is.na(protB)]

 If translate() already returns single letters, we can use them directly:
aaA <- protA
aaB <- protB

aa_levels <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*")
count_aa <- function(vec) as.numeric(table(factor(vec, levels = aa_levels)))

aaA_ct <- count_aa(aaA); aaB_ct <- count_aa(aaB)

 Exclude stop "*" for frequency denominator
keep <- aa_levels != "*"
aa_df <- data.frame(
  AA = aa_levels[keep],
  Freq_E_coli = aaA_ct[keep] / sum(aaA_ct[keep]),
  Freq_Saprospirales = aaB_ct[keep] / sum(aaB_ct[keep])
)
print(aa_df)

 Barplot (amino acids) with colour
barplot(rbind(aa_df$Freq_E_coli, aa_df$Freq_Saprospirales),
        beside = TRUE,
        names.arg = aa_df$AA,
        las = 2, ylab = "Frequency",
        main = "Amino-acid frequencies (translated CDS)",
        col = c("skyblue", "lightpink"))
legend("topright", legend = c("E. coli","Saprospirales"), bty = "n",
       fill = c("skyblue", "lightpink"))
grid()

# 5. Codon-usage bias-
STEP 0 — constants & helpers

VALID_BASES <- c("A","C","G","T")
STOP_CODONS <- c("TAA","TAG","TGA")

 All 64 DNA codons (TTT..GGG) and the 61 sense codons
.all_codons <- as.vector(outer(outer(VALID_BASES, VALID_BASES, paste0), VALID_BASES, paste0))
.SENSE_CODONS <- setdiff(.all_codons, STOP_CODONS)

 0.1 Clean a character vector of triplets: keep only A/C/G/T triplets
.clean_triplets <- function(x) {
  x <- toupper(x)
  x[grepl("^[ACGT]{3}$", x)]
}

 0.2 Standard genetic code (RNA codons with U)
GENETIC_CODE <- c(
  UUU="F", UUC="F", UUA="L", UUG="L",
  CUU="L", CUC="L", CUA="L", CUG="L",
  AUU="I", AUC="I", AUA="I", AUG="M",
  GUU="V", GUC="V", GUA="V", GUG="V",
  UCU="S", UCC="S", UCA="S", UCG="S",
  CCU="P", CCC="P", CCA="P", CCG="P",
  ACU="T", ACC="T", ACA="T", ACG="T",
  GCU="A", GCC="A", GCA="A", GCG="A",
  UAU="Y", UAC="Y", UAA="*", UAG="*",
  CAU="H", CAC="H", CAA="Q", CAG="Q",
  AAU="N", AAC="N", AAA="K", AAG="K",
  GAU="D", GAC="D", GAA="E", GAG="E",
  UGU="C", UGC="C", UGA="*", UGG="W",
  CGU="R", CGC="R", CGA="R", CGG="R",
  AGU="S", AGC="S", AGA="R", AGG="R",
  GGU="G", GGC="G", GGA="G", GGG="G"
)

 0.3 DNA (with T) -> AA using GENETIC_CODE (returns NA for unknown/stop)
codon2aa <- function(codon) {
  if (is.na(codon) || !grepl("^[ACGT]{3}$", codon)) return(NA_character_)
  key <- chartr("T", "U", codon)               # T->U
  aa  <- GENETIC_CODE[[key]]
  if (is.null(aa) || aa == "*") NA_character_ else aa
}


 STEP 1 — count sense codons genome-wide
 (frame +1, CDS trimmed to multiple of 3)

get_codon_counts <- function(cds_list) {
  counts <- setNames(integer(length(.SENSE_CODONS)), .SENSE_CODONS)
  for (seqv in cds_list) {
    n  <- length(seqv)
    n3 <- n - (n %% 3)
    if (n3 < 3) next
    s   <- toupper(seqv[1:n3])
    idx <- seq.int(1, n3, by = 3)
    cod <- paste0(s[idx], s[idx+1], s[idx+2])

     keep only strict A/C/G/T triplets and drop stops
    cod <- .clean_triplets(cod)
    if (!length(cod)) next
    cod <- cod[!(cod %in% STOP_CODONS)]
    if (!length(cod)) next

     tabulate only over the fixed set of sense codons (prevents unknown names)
    tt <- table(factor(cod, levels = .SENSE_CODONS))
    counts <- counts + as.integer(tt)
  }
  counts
}


 STEP 2 — RSCU and a simple bias score

rscu <- function(counts_named) {
  stopifnot(identical(sort(names(counts_named)), sort(.SENSE_CODONS)))
   map codon -> AA (NA for any issue, though names should be clean now)
  aa_map <- vapply(names(counts_named), codon2aa, FUN.VALUE = character(1))
  res <- setNames(rep(NA_real_, length(counts_named)), names(counts_named))
  for (aa in unique(aa_map[!is.na(aa_map)])) {
    idx <- which(aa_map == aa)
    tot <- sum(counts_named[idx])
    if (tot > 0) {
      k <- length(idx)                    # number of synonymous codons
      res[idx] <- counts_named[idx] / (tot / k)
    }
  }
  res
}

bias_score <- function(counts_named, rscu_vals) {
  ok <- !is.na(rscu_vals)
  w  <- counts_named[ok]; r <- rscu_vals[ok]
  if (sum(w) == 0) return(NA_real_)
  sum(w * abs(r - 1)) / sum(w)
}


 STEP 3 — GC content at 3rd codon position

gc3 <- function(cds_list) {
  gc <- 0L; tot <- 0L
  for (seqv in cds_list) {
    n  <- length(seqv)
    n3 <- n - (n %% 3)
    if (n3 < 3) next
    third <- toupper(seqv[seq(3, n3, by = 3)])
     count only strict A/C/G/T
    third <- third[third %in% VALID_BASES]
    gc  <- gc  + sum(third %in% c("G","C"))
    tot <- tot + length(third)
  }
  if (tot == 0) NA_real_ else gc / tot
}


 STEP 4 — run it

codonA <- get_codon_counts(cdsA)
codonB <- get_codon_counts(cdsB)

 ensure same order of names
codonB <- codonB[names(codonA)]

rscuA <- rscu(codonA)
rscuB <- rscu(codonB)

biasA <- bias_score(codonA, rscuA)
biasB <- bias_score(codonB, rscuB)

gc3A  <- gc3(cdsA)
gc3B  <- gc3(cdsB)


 STEP 5 — tidy table for your report

AA_for_codon <- vapply(names(codonA), codon2aa, FUN.VALUE = character(1))
q5_table <- data.frame(
  Codon = names(codonA),
  AA    = AA_for_codon,
  Count_E_coli        = as.integer(codonA),
  Count_Saprospirales = as.integer(codonB),
  RSCU_E_coli         = as.numeric(rscuA),
  RSCU_Saprospirales  = as.numeric(rscuB),
  row.names = NULL
)
print(head(q5_table, 12))
 write.csv(q5_table, "Q5_codon_usage_RSCU.csv", row.names = FALSE)

q5_bias <- data.frame(
  Organism = c(orgA_name, orgB_name),
  Bias_mean_abs_RSCU = c(biasA, biasB),
  GC3 = c(gc3A, gc3B)
)
print(q5_bias)


 STEP 6 — coloured RSCU plot

ord <- order(q5_table$AA, q5_table$Codon)
mat_rscu <- rbind(q5_table$RSCU_E_coli[ord], q5_table$RSCU_Saprospirales[ord])
lab <- paste0(q5_table$Codon[ord], "(", q5_table$AA[ord], ")")

op <- par(mar = c(10, 4, 4, 2) + 0.1)
barplot(mat_rscu, beside = TRUE, names.arg = lab, las = 2,
        ylab = "RSCU", main = "Codon usage bias (RSCU)",
        col = c("skyblue", "lightpink"))
abline(h = 1, lty = 2)
legend("topright", legend = c(orgA_name, orgB_name),
       bty = "n", fill = c("skyblue", "lightpink"))
grid()
par(op)

# 6. k-mer over/under-representation-
rgA_nam# Build AA character vectors (exclude NA and stop for counting)
protA <- unlist(lapply(cdsA, safe_translate), use.names = FALSE); protA <- protA[!is.na(protA)]
protB <- unlist(lapply(cdsB, safe_translate), use.names = FALSE); protB <- protB[!is.na(protB)]

collapse_aa <- function(prot_vec) unlist(strsplit(paste0(prot_vec, collapse = ""), split = "", fixed = TRUE), use.names = FALSE)
aaA <- collapse_aa(protA)
aaB <- collapse_aa(protB)

 Count kmers and compute expected under 0th-order AA model
kmers_counts <- function(aa_charvec, k) {
  aa <- aa_charvec[aa_charvec != "*"]
  n <- length(aa); if (n < k) return(integer(0))
  keys <- sapply(seq_len(n - k + 1), function(i) paste0(aa[i:(i+k-1)], collapse = ""))
  table(keys)
}
expected_counts <- function(aa_charvec, k) {
  aa <- aa_charvec[aa_charvec != "*"]
  n <- length(aa); if (n < k) return(list(obs=integer(0), exp=numeric(0)))
  freqs <- prop.table(table(aa))
  obs <- kmers_counts(aa_charvec, k)
  expv <- numeric(length(obs)); names(expv) <- names(obs)
  for (j in seq_along(obs)) {
    mer <- unlist(strsplit(names(obs)[j], ""))
    p <- prod(freqs[mer], na.rm = TRUE)
    expv[j] <- p * (n - k + 1)
  }
  list(obs = obs, exp = expv)
}

analyze_k <- function(aaA_char, aaB_char, k, ref_label="E_coli", test_label="Saprospirales") {
  aA <- expected_counts(aaA_char, k)
  aB <- expected_counts(aaB_char, k)
  all_k <- sort(unique(c(names(aA$obs), names(aB$obs))))
  obsA <- as.numeric(aA$obs[all_k]); obsA[is.na(obsA)] <- 0
  obsB <- as.numeric(aB$obs[all_k]); obsB[is.na(obsB)] <- 0
  expA <- as.numeric(aA$exp[all_k]); expA[is.na(expA)] <- 0
  expB <- as.numeric(aB$exp[all_k]); expB[is.na(expB)] <- 0
  eps <- 1e-9
  enrA <- (obsA + eps) / (expA + eps)
  enrB <- (obsB + eps) / (expB + eps)

  df <- data.frame(kmer = all_k, obs_A = obsA, exp_A = expA, enrich_A = enrA,
                               obs_B = obsB, exp_B = expB, enrich_B = enrB)
   Top 10 most over- and under-represented in Saprospirales
  ord_over  <- order(df$enrich_B, decreasing = TRUE)
  ord_under <- order(df$enrich_B, decreasing = FALSE)
  top_over  <- df[head(ord_over, 10), c("kmer","enrich_A","enrich_B","obs_A","obs_B")]
  top_under <- df[head(ord_under, 10), c("kmer","enrich_A","enrich_B","obs_A","obs_B")]

  cat("\n===== k =", k, " TOP 10 OVER in", test_label, "=====\n"); print(top_over, row.names = FALSE)
  cat("\n===== k =", k, " TOP 10 UNDER in", test_label, "=====\n"); print(top_under, row.names = FALSE)

   Scatter plot: enrichment comparison
  plot(df$enrich_A, df$enrich_B, pch = 16, cex = 0.6,
       xlab = paste0("Enrichment (", ref_label, ")"),
       ylab = paste0("Enrichment (", test_label, ")"),
       main = paste0("Protein k-mer enrichment (k=", k, ")"))
  abline(0,1,lty=2); grid()

   Barplots for tops (compare enrichments in both orgs)
  par(mfrow = c(1,2))
  barplot(rbind(top_over$enrich_A, top_over$enrich_B), beside = TRUE, las = 2,
          names.arg = top_over$kmer, ylab = "obs/exp",
          main = paste0("Top 10 OVER in ", test_label, " (k=", k, ")"))
  legend("topright", legend = c(ref_label, test_label), bty = "n", fill = grey.colors(2))
  barplot(rbind(top_under$enrich_A, top_under$enrich_B), beside = TRUE, las = 2,
          names.arg = top_under$kmer, ylab = "obs/exp",
          main = paste0("Top 10 UNDER in ", test_label, " (k=", k, ")"))
  legend("topright", legend = c(ref_label, test_label), bty = "n", fill = grey.colors(2))
  par(mfrow = c(1,1))

  invisible(list(all = df, top_over = top_over, top_under = top_under))
}

res_k3 <- analyze_k(aaA, aaB, k = 3, ref_label = orgA_name, test_label = orgB_name)
res_k4 <- analyze_k(aaA, aaB, k = 4, ref_label = orgA_name, test_label = orgB_name)
res_k5 <- analyze_k(aaA, aaB, k = 5, ref_label = orgA_name, test_label = orgB_name)


# Purpose of each script
These are the part2_download_sequences. The R tool gets the CDS and protein FASTA files for both the species.  These are the part2_sequence_summaries. The R script counts the genes, figures out the length of all the coding DNA, and makes boxplots that show how the CDS lengths are distributed.  The frequency and codon for 05_part2_comp. The R script figures out the rates of nucleotides and amino acids, shows them with barplots, and makes codon usage tables with bias measures like RSCU values.  Last but not least, 06_part2_kmer_enrichment. R finds over- and under-represented protein k-mers (length 3–5) in the Saprospirales bacteria, compares them to E. coli, and makes barplots that show changes in sequence bias.

# Troubleshooting
If a series of downloads fails, delete any that aren't finished. file and run 03_part2_download_sequences.R again.  Put in TinyTeX or knit to HTML format instead if the report won't knit because of LaTeX problems.  Before running the scripts, make sure that all of the relevant packages have been loaded.  If memory errors happen while FASTA parsing, close all other programs and try again; the sequences are pretty small and usually load without any problems.

# Data Citation
The CDS data that was used in this study came from NCBI and Ensembl Bacteria.  All of the datasets are free to use and are meant to be used for learning.  Use the phrase "Ensembl Bacteria release 62 and NCBI RefSeq data for E. coli and Saprospirales bacterium" when you cite them.

# Academic Integrity
This library holds the results of independent research. There is no copied code or work, and some AI tools are used. To keep academic ethics and make sure that the records can be used again, they have not been changed. Structure and comments are included to help others understand the analysis process in a clear way.

# License and Authorship
Author- Suhan Mhatre
License- MIT
