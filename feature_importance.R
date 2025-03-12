# Identify genes that best distinguish specific cancer types from all others
# This script extracts feature importance from the Random Forest model

# Load required libraries
library(dplyr)
library(ggplot2)
library(randomForest)
library(caret)
library(tibble)  # Add this library for column_to_rownames function
library(pheatmap)
library(tidyr)  # Explicitly load for pivot_longer

# Load the pre-trained Random Forest model
rf_model <- readRDS("random_forest_model.rds")

# Load processed data
processed_data <- readRDS("ml_processed_data.rds")
test_features <- processed_data$test_features
test_target <- processed_data$test_target

# Get variable importance from the model
var_importance <- varImp(rf_model)
importance_df <- var_importance$importance
importance_df$Gene <- rownames(importance_df)

# Rearrange the data for easier handling
importance_long <- pivot_longer(
  importance_df, 
  cols = -Gene, 
  names_to = "Cancer_Type", 
  values_to = "Importance"
)

# Sort and get the top genes for each cancer type
top_genes <- importance_long %>%
  group_by(Cancer_Type) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Save the top genes for each cancer type
write.csv(top_genes, "top_genes_by_cancer_type.csv", row.names = FALSE)

# Get overall top genes across all cancer types
overall_importance <- importance_df %>%
  mutate(Overall_Importance = rowMeans(select(., -Gene))) %>%
  arrange(desc(Overall_Importance))

# Save overall top genes
write.csv(overall_importance[1:50, ], "top_50_genes_overall.csv", row.names = FALSE)

# Visualize top 20 genes overall
top_20_overall <- overall_importance[1:20, ]
png("top_20_genes_overall.png", width = 1000, height = 600)
ggplot(top_20_overall, aes(x = reorder(Gene, Overall_Importance), y = Overall_Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Genes with Highest Overall Importance",
       x = "Gene", y = "Average Importance") +
  theme_minimal()
dev.off()

# Select one specific cancer type for detailed analysis
# Choose the cancer type with highest prevalence or of particular interest
cancer_types <- unique(importance_long$Cancer_Type)
specific_cancer <- cancer_types[1]  # Replace with specific cancer of interest if needed

# Get top genes for the chosen cancer type
specific_cancer_genes <- importance_long %>%
  filter(Cancer_Type == specific_cancer) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = 20)

# Visualize top genes for the specific cancer type
png(paste0("top_genes_for_", specific_cancer, ".png"), width = 1000, height = 600)
ggplot(specific_cancer_genes, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "darkred") +
  coord_flip() +
  labs(title = paste("Top 20 Genes for Distinguishing", specific_cancer),
       x = "Gene", y = "Importance") +
  theme_minimal()
dev.off()

# Create a heatmap of top 30 genes across all cancer types
top_30_genes <- overall_importance$Gene[1:30]

# Alternative approach without column_to_rownames
importance_heatmap_data <- as.matrix(importance_df[importance_df$Gene %in% top_30_genes, -which(names(importance_df) == "Gene")])
rownames(importance_heatmap_data) <- importance_df$Gene[importance_df$Gene %in% top_30_genes]

# Order genes by overall importance
gene_means <- rowMeans(importance_heatmap_data)
importance_heatmap_data <- importance_heatmap_data[order(gene_means, decreasing = TRUE), ]

# Create heatmap using pheatmap
png("gene_importance_heatmap.png", width = 1200, height = 800, res = 120)
pheatmap(importance_heatmap_data,
         main = "Importance of Top 30 Genes Across Cancer Types",
         color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         fontsize_row = 8,
         fontsize_col = 8)
dev.off()

cat("Analysis complete. Files saved:\n",
    "- top_genes_by_cancer_type.csv: Top 10 genes for each cancer type\n",
    "- top_50_genes_overall.csv: Top 50 genes overall\n",
    "- top_20_genes_overall.png: Visualization of top 20 genes\n",
    "- top_genes_for_[cancer_type].png: Visualization for specific cancer\n",
    "- gene_importance_heatmap.png: Heatmap of top 30 genes\n")