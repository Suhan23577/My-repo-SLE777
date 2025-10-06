# My-repo-SLE777
ASSIGNMENT 4 (Part 1)- Data Wrangling, Analysis and Plotting

# Overview
This repository has R tools and a R Markdown report that can be used with RStudio and GitHub to organize data, display it, and analyze sequence variation.  Parts 1 and 2 are both covered. Part 1 has points 1 through 10; Part 2 has points 1 through 6.  In Part 1, the tools bring in and look at two sets of data (gene_expression.tsv and growth_data.csv), make plots, create descriptive statistics, and test hypotheses.  In Part 2, they get coding DNA sequences (CDS) for E. coli and Saprospirales bacteria and then figure out data at the sequence level, like gene counts, nucleotide rates, codon usage bias, and k-mer enrichment.  All of the code works perfectly, produces results that can be repeated, and can be put together into a single, complete report using R Markdown.

# Repository Structure
The repository is set up so that it is easy to understand and use again.  Files like gene_expression.tsv and growth_data.csv are stored in the data/ folder.  There are two subfolders, part1/ and part2/, inside the outputs/ folder that hold all the figures and tables that were made.  There are several R scripts in the scripts/ directory. Each one does a different part of the project, from importing and organizing data to visualizing it and analyzing sequences.  These tools are put together in the main R Markdown file, report.Rmd, to make the end report.  For faster work in RStudio and GitHub, supporting files like.Rproj and.gitignore are included.

# Data Sources
download and save files locally:
download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/main/gene_expression.tsv",
              destfile = "gene_expression.tsv")

download.file("https://raw.githubusercontent.com/ghazkha/Assessment4/main/growth_data.csv",
              destfile = "growth_data.csv")

now you can read them:
gene_data <- read.table("gene_expression.tsv", header = TRUE, sep = "\t", row.names = 1)
growth_data <- read.csv("growth_data.csv")

# How to Run
1. Read file with gene IDs as row names; show first six genes-

Import the tab-separated file
gene_data <- read.table("gene_expression.tsv", 
                        header = TRUE, 
                        sep = "\t", 
                        row.names = 1)

View the first 6 genes (rows)
head(gene_data)

2. Add a mean column; show first six-

Add a column for the mean expression value across all samples
gene_data$mean_expression <- rowMeans(gene_data)

Display the first six genes again
head(gene_data)

3. List top 10 genes by mean expression-

Sort by mean expression and display the top 10
top10 <- head(gene_data[order(-gene_data$mean_expression), ], 10)
top10

4. Count genes with mean < 10-

Count how many genes have mean expression below 10
low_genes <- sum(gene_data$mean_expression < 10)
low_genes

5. Histogram of mean values-

Plot histogram of mean expression values
hist(gene_data$mean_expression,
     main = "Distribution of Mean Gene Expression",
     xlab = "Mean Expression Value",
     ylab = "Frequency",
     col = "lavender",
     border = "pink")

6. Import CSV; print column names-

Read the CSV data
growth_data <- read.csv("growth_data.csv", header = TRUE)

Display column names
colnames(growth_data)

7. Mean & SD at start (2005) and end (2020) by site-

aggregate(cbind(Circumf_2005_cm, Circumf_2020_cm) ~ Site, 
          data = growth_data,
          FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                              sd = sd(x, na.rm = TRUE)))
Make sure columns are numeric (safe if they already are)
growth_data$Circumf_2005_cm <- as.numeric(growth_data$Circumf_2005_cm)
growth_data$Circumf_2020_cm <- as.numeric(growth_data$Circumf_2020_cm)

Compute mean & sd by Site
agg <- aggregate(cbind(Circumf_2005_cm, Circumf_2020_cm) ~ Site,
                 data = growth_data,
                 FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                     sd   = sd(x,   na.rm = TRUE)))

Unnest the matrix columns into a clean data frame
out <- data.frame(
  Site       = agg$Site,
  Start_mean = agg$Circumf_2005_cm[, "mean"],
  Start_sd   = agg$Circumf_2005_cm[, "sd"],
  End_mean   = agg$Circumf_2020_cm[, "mean"],
  End_sd     = agg$Circumf_2020_cm[, "sd"],
  row.names = NULL
)

FORCE the display (works in scripts, Rmd, and notebooks)
print(out)

8. Boxplots at start vs end by site-

boxplot(Circumf_2005_cm ~ Site, data = growth_data,
        main = "Tree Circumference in 2005 (Start)",
        ylab = "Circumference (cm)", col = "pink")

boxplot(Circumf_2020_cm ~ Site, data = growth_data,
        main = "Tree Circumference in 2020 (End)",
        ylab = "Circumference (cm)", col = "lavender")

9. Mean growth over last 10 years at each site-

growth_data$growth_2010_2020 <- growth_data$Circumf_2020_cm - growth_data$Circumf_2010_cm
aggregate(growth_2010_2020 ~ Site, data = growth_data, mean)

10. t-test: is 10-year growth different between sites-

t.test(growth_2010_2020 ~ Site, data = growth_data)


# Purpose of each script
The 00_setup.R script gets the environment ready by making the files it needs, checking for packages, and setting up functions that will be useful.  This script is called 01_part1_gene_expression. R reads the gene expression file, writes the gene identifiers as row names, finds the ten genes with the highest mean expression values, counts the genes whose mean expression value is less than 10, and makes a plot of those mean expression values.  The code for 02_part1_growth_analysis. The tree growth data is imported into R, column names are listed, and the mean and standard deviation for circumference are found at both sites (the beginning and end of the study). Boxplots are made to show how growth has changed at each site, the mean 10-year growth is found, and a t-test is used to compare growth between the control and treatment sites.








