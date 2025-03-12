# TCGA Cancer Type Prediction Model
# This script builds and evaluates machine learning models to predict cancer type from gene expression

# Load required libraries
library(dplyr)
library(caret)
library(randomForest)
library(e1071)
library(ggplot2)
library(doParallel)

# Optional: Register parallel processing to speed up training
registerDoParallel(cores = 4)  # Adjust based on your computer's capabilities

# Load the dataset
tcga_data <- readRDS("Dataset.rds")

# Data preparation
# Extract features (gene expression) and target (cancer type)
features <- tcga_data[, 3:ncol(tcga_data)]  # Assuming first two columns are Sample and Cancertype
target <- tcga_data$Cancertype
target <- as.factor(target)  # Ensure cancer type is treated as categorical

# Check class distribution
class_counts <- table(target)
print("Cancer type distribution:")
print(class_counts)

# Create train/test split (80/20)
set.seed(42)  # For reproducibility
train_index <- createDataPartition(target, p = 0.8, list = FALSE)
train_features <- features[train_index, ]
train_target <- target[train_index]
test_features <- features[-train_index, ]
test_target <- target[-train_index]

# Check dimensions of training and test sets
cat("Training set dimensions:", dim(train_features), "\n")
cat("Test set dimensions:", dim(test_features), "\n")

# Feature selection
# Optional: Select top features to reduce dimensionality and training time
if(ncol(features) > 500) {
  message("Performing feature selection...")
  
  # Calculate variance of each gene
  gene_var <- apply(train_features, 2, var)
  
  # Select top 500 most variable genes
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:500]
  
  # Subset data to include only top genes
  train_features <- train_features[, top_genes]
  test_features <- test_features[, top_genes]
  
  cat("Reduced features to top 500 most variable genes\n")
}

# Save the processed data for potential future use
processed_data <- list(
  train_features = train_features,
  train_target = train_target,
  test_features = test_features,
  test_target = test_target
)
saveRDS(processed_data, "ml_processed_data.rds")

# 1. Random Forest Model
message("Training Random Forest model...")
# Define training control - use cross-validation
train_control <- trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  classProbs = TRUE,
  savePredictions = "final"
)

# Train the Random Forest model
rf_model <- train(
  x = train_features,
  y = train_target,
  method = "rf",
  trControl = train_control,
  ntree = 100,  # Reduce for faster training
  importance = TRUE
)

# Save the trained model
saveRDS(rf_model, "random_forest_model.rds")

# Make predictions on test set
rf_predictions <- predict(rf_model, test_features)

# Evaluate Random Forest model
rf_confusion <- confusionMatrix(rf_predictions, test_target)
print("Random Forest Results:")
print(rf_confusion)

# Extract class-specific metrics
rf_by_class <- rf_confusion$byClass
rf_metrics <- data.frame(
  Cancer_Type = rownames(rf_by_class),
  Precision = rf_by_class[, "Pos Pred Value"],
  Recall = rf_by_class[, "Sensitivity"],
  F1_Score = rf_by_class[, "F1"],
  Specificity = rf_by_class[, "Specificity"]
)

# Save detailed metrics
write.csv(rf_metrics, "random_forest_metrics.csv", row.names = FALSE)

# Get variable importance
rf_importance <- varImp(rf_model)
top_genes <- rownames(rf_importance$importance)[1:20]

# 2. SVM Model
message("Training SVM model...")
# Train the SVM model
svm_model <- train(
  x = train_features,
  y = train_target,
  method = "svmLinear",
  trControl = train_control
)

# Save the trained model
saveRDS(svm_model, "svm_model.rds")

# Make predictions on test set
svm_predictions <- predict(svm_model, test_features)

# Evaluate SVM model
svm_confusion <- confusionMatrix(svm_predictions, test_target)
print("SVM Results:")
print(svm_confusion)

# Extract class-specific metrics
svm_by_class <- svm_confusion$byClass
svm_metrics <- data.frame(
  Cancer_Type = rownames(svm_by_class),
  Precision = svm_by_class[, "Pos Pred Value"],
  Recall = svm_by_class[, "Sensitivity"],
  F1_Score = svm_by_class[, "F1"],
  Specificity = svm_by_class[, "Specificity"]
)

# Save detailed metrics
write.csv(svm_metrics, "svm_metrics.csv", row.names = FALSE)

# Compare models
model_comparison <- data.frame(
  Model = c("Random Forest", "SVM"),
  Accuracy = c(rf_confusion$overall["Accuracy"], svm_confusion$overall["Accuracy"]),
  Kappa = c(rf_confusion$overall["Kappa"], svm_confusion$overall["Kappa"])
)

print("Model Comparison:")
print(model_comparison)
write.csv(model_comparison, "model_comparison.csv", row.names = FALSE)

# Create visualizations
# 1. Accuracy comparison
png("model_accuracy_comparison.png", width = 800, height = 600)
ggplot(model_comparison, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity") +
  labs(title = "Model Accuracy Comparison", y = "Accuracy") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1)) +
  geom_text(aes(label = sprintf("%.2f%%", Accuracy*100)), vjust = -0.5)
dev.off()

# 2. Class-specific metrics for Random Forest
long_rf_metrics <- rf_metrics %>%
  tidyr::pivot_longer(cols = c("Precision", "Recall", "F1_Score"), 
                      names_to = "Metric", values_to = "Value")

# Sort by F1 score for better visualization
top_cancer_types <- rf_metrics %>%
  arrange(desc(F1_Score)) %>%
  head(15) %>%
  pull(Cancer_Type)

filtered_metrics <- long_rf_metrics %>%
  filter(Cancer_Type %in% top_cancer_types)

png("rf_metrics_by_cancer_type.png", width = 1200, height = 800)
ggplot(filtered_metrics, aes(x = reorder(Cancer_Type, Value), y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Random Forest Performance by Cancer Type",
       x = "Cancer Type", y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 1))
dev.off()

# 3. Confusion matrix heatmap for RF (simplified for visibility - top 10 cancer types)
top_10_cancer_types <- names(sort(table(test_target), decreasing = TRUE))[1:10]
top_10_indices <- which(rownames(rf_confusion$table) %in% top_10_cancer_types)
conf_matrix_subset <- rf_confusion$table[top_10_indices, top_10_indices]

# Normalize by row to show prediction patterns
conf_matrix_norm <- sweep(conf_matrix_subset, 1, rowSums(conf_matrix_subset), "/")

# Create data frame for ggplot
conf_df <- as.data.frame(as.table(conf_matrix_norm))
names(conf_df) <- c("Actual", "Predicted", "Proportion")

png("rf_confusion_matrix_heatmap.png", width = 1000, height = 900)
ggplot(conf_df, aes(x = Predicted, y = Actual, fill = Proportion)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Random Forest Confusion Matrix (Top 10 Cancer Types)",
       x = "Predicted", y = "Actual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Print summary message
cat("\nMachine learning analysis complete.\n",
    "Random Forest accuracy: ", round(rf_confusion$overall["Accuracy"]*100, 2), "%\n",
    "SVM accuracy: ", round(svm_confusion$overall["Accuracy"]*100, 2), "%\n",
    "Results saved to CSV files and visualizations.\n")