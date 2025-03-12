# TCGA Dataset Analysis
# This script analyzes TCGA gene expression data to:
# 1. Count samples per cancer type
# 2. Calculate mean and standard deviation for each gene
# 3. Identify the top 5 most variable genes

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the dataset
tcga_data <- readRDS("Dataset.rds")

# 1. Number of samples per cancer type
cancer_counts <- table(tcga_data$Cancertype)
print("Number of samples per cancer type:")
print(cancer_counts)

# Create a barplot of sample counts by cancer type
png("cancer_type_counts.png", width=800, height=600)
barplot(cancer_counts, 
        main="Number of Samples by Cancer Type", 
        xlab="Cancer Type", 
        ylab="Number of Samples",
        las=2,  # Rotate labels
        cex.names=0.7)  # Adjust label size
dev.off()

# 2. Mean and standard deviation of expression levels for each gene
# Exclude non-gene columns (Sample and Cancertype)
gene_columns <- 3:ncol(tcga_data)  # Assuming first two columns are Sample and Cancertype
gene_names <- names(tcga_data)[gene_columns]

gene_stats <- data.frame(
  gene = gene_names,
  mean = sapply(tcga_data[, gene_columns], mean),
  sd = sapply(tcga_data[, gene_columns], sd)
)

# Save gene statistics to a CSV file
write.csv(gene_stats, "gene_statistics.csv", row.names = FALSE)
print("First few rows of gene statistics:")
print(head(gene_stats))

# 3. Top 5 most variable genes across all samples
gene_variance <- sapply(tcga_data[, gene_columns], var)
top_variable_genes <- names(sort(gene_variance, decreasing = TRUE))[1:5]
print("Top 5 most variable genes:")
print(top_variable_genes)

# Show the variance values for these top genes
top_gene_variances <- sort(gene_variance, decreasing = TRUE)[1:5]
variance_df <- data.frame(Gene = names(top_gene_variances), Variance = top_gene_variances)
print(variance_df)

# Save top variable genes info to a CSV file
write.csv(variance_df, "top_variable_genes.csv", row.names = FALSE)

# Create a bar plot of the top variable genes
png("top_variable_genes.png", width=600, height=400)
barplot(top_gene_variances, 
        names.arg=names(top_gene_variances),
        main="Top 5 Most Variable Genes",
        xlab="Gene",
        ylab="Variance",
        col="steelblue")
dev.off()

# Summary of results
cat("\nAnalysis Summary:\n")
cat("1. Found", length(cancer_counts), "cancer types with a total of", nrow(tcga_data), "samples\n")
cat("2. Calculated statistics for", nrow(gene_stats), "genes\n")
cat("3. Identified top variable genes:", paste(top_variable_genes, collapse=", "), "\n")