# TCGA Clustering Analysis
# This script performs hierarchical and k-means clustering on TCGA data and visualizes the results

# Load required libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Rtsne)
library(factoextra)
library(RColorBrewer)

# Load the dataset
tcga_data <- readRDS("Dataset.rds")

# Extract gene expression data (excluding Sample and Cancertype columns)
gene_columns <- 3:ncol(tcga_data)
expression_matrix <- as.matrix(tcga_data[, gene_columns])
rownames(expression_matrix) <- rownames(tcga_data)

# Create a color vector for cancer types that works with more than 12 types
# First check how many cancer types we have
n_cancer_types <- length(unique(tcga_data$Cancertype))
cancer_types <- unique(tcga_data$Cancertype)

# Generate colors - using a different approach to handle more than 12 types
if (n_cancer_types <= 12) {
  colors <- brewer.pal(n_cancer_types, "Set3")
} else {
  # For more than 12 types, we can combine multiple palettes or use a different approach
  colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_cancer_types)
}

# Create the named vector
cancer_colors <- setNames(colors, cancer_types)

# 1. Hierarchical Clustering
# Compute distance matrix (using correlation-based distance for gene expression)
message("Computing distance matrix...")
dist_matrix <- dist(expression_matrix, method = "euclidean")

# Perform hierarchical clustering
message("Performing hierarchical clustering...")
hc <- hclust(dist_matrix, method = "ward.D2")

# Create a cluster assignment with a reasonable number of clusters
k <- length(unique(tcga_data$Cancertype))  # Use number of cancer types as k
cluster_assignments <- cutree(hc, k = k)

# Create data frame with cluster assignments and cancer types
cluster_data <- data.frame(
  Sample = rownames(expression_matrix),
  Cancer_Type = tcga_data$Cancertype,
  Cluster = as.factor(cluster_assignments)
)

# Save cluster assignments
write.csv(cluster_data, "hierarchical_cluster_assignments.csv", row.names = FALSE)

# Visualize dendrogram with colored cancer types
png("hierarchical_clustering_dendrogram.png", width = 1500, height = 1000, res = 120)
plot(hc, labels = FALSE, main = "Hierarchical Clustering Dendrogram", sub = "", xlab = "", ylab = "Height")
rect.hclust(hc, k = k, border = "red")
# Add a simple legend (with fewer items if too many cancer types)
if (n_cancer_types <= 20) {
  legend("topright", legend = cancer_types, 
         fill = cancer_colors, 
         cex = 0.7, title = "Cancer Type")
}
dev.off()

# Create a confusion matrix to compare clusters vs cancer types
confusion_matrix <- table(cluster_data$Cluster, cluster_data$Cancer_Type)
write.csv(confusion_matrix, "cluster_vs_cancer_confusion_matrix.csv")

# 2. Dimensionality reduction for visualization
# PCA with handling for zero-variance genes
message("Performing PCA...")
# Identify and remove constant/zero columns before PCA
gene_variances <- apply(expression_matrix, 2, var)
zero_var_genes <- which(gene_variances == 0 | is.na(gene_variances))

if(length(zero_var_genes) > 0) {
  message(paste("Removing", length(zero_var_genes), "genes with zero variance before PCA"))
  filtered_expression_matrix <- expression_matrix[, -zero_var_genes]
} else {
  filtered_expression_matrix <- expression_matrix
}

# Now perform PCA on the filtered matrix
pca_result <- prcomp(filtered_expression_matrix, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$Cancer_Type <- tcga_data$Cancertype

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
var_explained_labels <- paste0("PC", 1:2, " (", round(var_explained[1:2] * 100, 1), "%)")

# PCA plot
png("pca_plot.png", width = 1200, height = 1000, res = 120)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cancer_Type)) +
  geom_point(alpha = 0.7) +
  labs(x = var_explained_labels[1], y = var_explained_labels[2],
       title = "PCA of TCGA Gene Expression Data", 
       subtitle = "Colored by cancer type") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(values = cancer_colors)
dev.off()

# t-SNE (may take some time for large datasets)
message("Performing t-SNE...")
set.seed(42)  # For reproducibility
# Use the filtered matrix for t-SNE as well
tsne_result <- Rtsne(filtered_expression_matrix, dims = 2, 
                    perplexity = min(30, nrow(filtered_expression_matrix)/5), 
                    check_duplicates = FALSE, pca = TRUE, 
                    pca_center = TRUE, pca_scale = TRUE,
                    max_iter = 1000)
tsne_data <- as.data.frame(tsne_result$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")
tsne_data$Cancer_Type <- tcga_data$Cancertype

# t-SNE plot
png("tsne_plot.png", width = 1200, height = 1000, res = 120)
ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = Cancer_Type)) +
  geom_point(alpha = 0.7) +
  labs(title = "t-SNE of TCGA Gene Expression Data",
       subtitle = "Colored by cancer type") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_manual(values = cancer_colors)
dev.off()

# 3. k-means clustering (use filtered expression matrix here too)
message("Performing k-means clustering...")
set.seed(42)
k_means_result <- kmeans(filtered_expression_matrix, centers = k, nstart = 25)

# Add k-means cluster to data
kmeans_data <- data.frame(
  Sample = rownames(filtered_expression_matrix),
  Cancer_Type = tcga_data$Cancertype,
  Cluster = as.factor(k_means_result$cluster)
)

# Save k-means cluster assignments
write.csv(kmeans_data, "kmeans_cluster_assignments.csv", row.names = FALSE)

# K-means confusion matrix
kmeans_confusion <- table(kmeans_data$Cluster, kmeans_data$Cancer_Type)
write.csv(kmeans_confusion, "kmeans_vs_cancer_confusion_matrix.csv")

# Visualize k-means with PCA
pca_with_kmeans <- pca_data
pca_with_kmeans$Cluster <- as.factor(k_means_result$cluster)

png("pca_with_kmeans.png", width = 1200, height = 1000, res = 120)
ggplot(pca_with_kmeans, aes(x = PC1, y = PC2, color = Cluster, shape = Cancer_Type)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = var_explained_labels[1], y = var_explained_labels[2],
       title = "PCA with K-means Clustering",
       subtitle = "Shape represents cancer type, color represents cluster") +
  theme_minimal() +
  guides(shape = guide_legend(override.aes = list(size = 3)))
dev.off()

# 4. Correlation between cancer types
# Calculate average expression per cancer type
# We can use the original expression matrix here
avg_expression_by_cancer <- matrix(0, nrow = length(cancer_types), ncol = ncol(expression_matrix))
rownames(avg_expression_by_cancer) <- cancer_types
colnames(avg_expression_by_cancer) <- colnames(expression_matrix)

for(i in seq_along(cancer_types)) {
  cancer_samples <- tcga_data$Cancertype == cancer_types[i]
  avg_expression_by_cancer[i,] <- colMeans(expression_matrix[cancer_samples,])
}

# Calculate correlation between cancer types
cancer_correlation <- cor(t(avg_expression_by_cancer))

# Heatmap of cancer type correlations
png("cancer_correlation_heatmap.png", width = 1000, height = 900, res = 120)
pheatmap(cancer_correlation, 
         main = "Correlation Between Cancer Types Based on Gene Expression",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = (n_cancer_types <= 15),  # Only show numbers if not too many types
         number_format = "%.2f",
         fontsize_number = 8)
dev.off()

# 5. Analysis summary
# Calculate agreement between clustering and cancer types
# For hierarchical clustering
hierarchical_agreement <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
# For k-means
kmeans_agreement <- sum(diag(kmeans_confusion)) / sum(kmeans_confusion)

# Create summary file
cat("Clustering Analysis Summary\n\n",
    "Number of samples: ", nrow(tcga_data), "\n",
    "Number of genes: ", ncol(expression_matrix), "\n",
    "Number of cancer types: ", length(cancer_types), "\n\n",
    "Hierarchical Clustering:\n",
    "- Agreement with cancer types: ", round(hierarchical_agreement * 100, 2), "%\n\n",
    "K-means Clustering:\n",
    "- Agreement with cancer types: ", round(kmeans_agreement * 100, 2), "%\n\n",
    "Interpretation:\n",
    "Based on the clustering results, we can observe whether cancers with similar tissue origins cluster together.\n",
    "The confusion matrices show how samples from each cancer type are distributed across clusters.\n",
    "The correlation heatmap shows similarities between different cancer types based on gene expression profiles.\n",
    file = "clustering_analysis_summary.txt")

message("Clustering analysis complete. Results saved to files.")