# TCGA Differential Expression Analysis (Complete)
# This script identifies differentially expressed genes between cancer types and creates visualizations

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)  # For heatmap visualization
library(RColorBrewer)  # For color palettes

# Load the dataset (if not already loaded)
tcga_data <- readRDS("Dataset.rds")

# Select 5 cancer types for comparison
# First check how many samples we have per cancer type
cancer_counts <- table(tcga_data$Cancertype)
print("Available cancer types:")
print(cancer_counts)

# Select the 5 cancer types with most samples (you can adjust this if needed)
top_5_cancers <- names(sort(cancer_counts, decreasing = TRUE)[1:5])
print("Selected cancer types for comparison:")
print(top_5_cancers)

# Filter data to include only these 5 cancer types
selected_data <- tcga_data[tcga_data$Cancertype %in% top_5_cancers, ]

# Perform differential expression analysis
gene_columns <- 3:ncol(tcga_data)
gene_names <- names(tcga_data)[gene_columns]
p_values <- numeric(length(gene_columns))

# Create a matrix to store mean expression values for each cancer type
mean_expression <- matrix(0, nrow = length(gene_columns), ncol = length(top_5_cancers))
colnames(mean_expression) <- top_5_cancers

for(i in seq_along(gene_columns)) {
  # Kruskal-Wallis test for each gene
  test_result <- kruskal.test(selected_data[[gene_columns[i]]] ~ selected_data$Cancertype)
  p_values[i] <- test_result$p.value
  
  # Calculate mean expression in each cancer type
  for(j in seq_along(top_5_cancers)) {
    cancer_type <- top_5_cancers[j]
    cancer_samples <- selected_data$Cancertype == cancer_type
    mean_expression[i, j] <- mean(selected_data[cancer_samples, gene_columns[i]])
  }
}

# Adjust p-values for multiple testing
adjusted_p_values <- p.adjust(p_values, method = "BH")

# Calculate overall variance for each gene (to measure effect size)
gene_variance <- apply(mean_expression, 1, var)
fold_change_max <- apply(mean_expression, 1, function(x) max(x) / (min(x) + 0.01))  # Adding 0.01 to avoid division by zero

# Create results dataframe with mean expression values for each cancer type
diff_exp_results <- data.frame(
  gene = gene_names,
  p_value = p_values,
  adjusted_p_value = adjusted_p_values,
  variance_between_cancers = gene_variance,
  max_fold_change = fold_change_max
)

# Add mean expression for each cancer type
for(j in seq_along(top_5_cancers)) {
  cancer_type <- top_5_cancers[j]
  diff_exp_results[[paste0("mean_", cancer_type)]] <- mean_expression[, j]
}

# Sort by adjusted p-value first, then by variance
diff_exp_results <- diff_exp_results[order(diff_exp_results$adjusted_p_value, -diff_exp_results$variance_between_cancers), ]

# Get top 10 differentially expressed genes
top_10_genes <- diff_exp_results$gene[1:10]
print("Top 10 differentially expressed genes:")
print(head(diff_exp_results[1:10, ], 10))

# Save complete results to CSV
write.csv(diff_exp_results, "differential_expression_results_complete.csv", row.names = FALSE)

# Create a heatmap of the top 10 genes
# Extract expression data for top genes and convert to matrix
expression_data <- as.matrix(selected_data[, top_10_genes])
rownames(expression_data) <- rownames(selected_data)

# Create annotation data frame for heatmap
annotation_df <- data.frame(
  CancerType = selected_data$Cancertype
)
rownames(annotation_df) <- rownames(selected_data)

# Create color palette for cancer types
cancer_colors <- list(
  CancerType = setNames(
    brewer.pal(length(top_5_cancers), "Set1"),
    top_5_cancers
  )
)

# Generate heatmap
png("top_genes_heatmap.png", width = 1200, height = 800, res = 120)
pheatmap(
  expression_data,
  annotation_row = annotation_df,
  annotation_colors = cancer_colors,
  main = "Expression of Top 10 Differentially Expressed Genes",
  fontsize = 8,
  fontsize_row = 6,
  show_rownames = FALSE,  # Too many rows to show names
  scale = "row",  # Scale by row (z-score)
  clustering_method = "ward.D2"
)
dev.off()

# Create a boxplot for each of the top 10 genes
dir.create("gene_boxplots", showWarnings = FALSE)

for(gene in top_10_genes) {
  plot_data <- data.frame(
    Cancer_Type = selected_data$Cancertype,
    Expression = selected_data[[gene]]
  )
  
  p <- ggplot(plot_data, aes(x = Cancer_Type, y = Expression, fill = Cancer_Type)) +
    geom_boxplot() +
    labs(title = paste("Expression of", gene, "Across Cancer Types"),
         x = "Cancer Type", y = "Expression Level") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(paste0("gene_boxplots/", gene, "_boxplot.png"), p, width = 8, height = 6)
}

# Create a summary heatmap with average expression by cancer type for top 10 genes
# We already have mean_expression matrix, but let's extract just the top 10 genes
top_10_indices <- match(top_10_genes, gene_names)
avg_expression <- t(mean_expression[top_10_indices, ])  # Transpose to have genes as columns

# Generate summary heatmap
png("avg_expression_heatmap.png", width = 800, height = 600, res = 100)
pheatmap(
  avg_expression,
  main = "Average Expression of Top 10 Genes by Cancer Type",
  fontsize = 10,
  cluster_rows = FALSE,
  cluster_cols = TRUE
)
dev.off()

# Create a correlation heatmap of the top 10 genes
gene_correlation <- cor(expression_data)
png("gene_correlation_heatmap.png", width = 800, height = 800, res = 100)
pheatmap(
  gene_correlation,
  main = "Correlation Between Top 10 Differentially Expressed Genes",
  fontsize = 10
)
dev.off()

print("Analysis complete. Output files saved.")