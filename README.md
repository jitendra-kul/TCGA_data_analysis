# TCGA_data_analysis

Questions:
Load the TCGA dataset and provide a summary including:
   1. The number of samples per cancer type
   2. The mean and standard deviation of expression levels for each gene
   3. The top 5 most variable genes across all samples

Another question is"Identify the top 10 differentially expressed genes between 5 cancer types using a Kruskal-Wallis test test. Visualize the expression of these genes with a heatmap. (you can make multiple visualization plots)"

![avg_expression_heatmap](https://github.com/user-attachments/assets/b13e52e4-c403-4518-a3c0-bd359a47a2ae)
![top_genes_heatmap](https://github.com/user-attachments/assets/3c52c2a1-6870-414a-84c1-43013d4b396b)


Perform hierarchical clustering or k-means clustering on the dataset using all genes.
   
•	Visualize the clustering with a dendrogram or t-SNE/PCA plot.

•	Interpret whether cancers with similar tissue origins cluster together


Here I performed clustering to determine if cancers with similar tissue origins have similar gene expression patterns.

Overview of the dataset: From the summary file which will be generated after a user will run the task3 code: it shows 9263 samples, 716 genes and 24 cancer types.

![Task_3_cancer_correlation_heatmap](https://github.com/user-attachments/assets/0f5513e1-2b48-4695-b7bc-9e6c60bd298e)

This heatmap shows correlations between different cancer types based on their gene expression profiles. Red indicates high correlation (similarity) while blue indicates low correlation. The hierarchical clustering (dendrogram) at the top groups cancer types with similar expression patterns. Note how cancers from similar tissue origins cluster together: kidney cancers (KIRC, KIRP, KICH), brain tumors (GBM, LGG), and squamous cell carcinomas (BLCA, LUSC, HNSC, CESC) form distinct groups. This provides strong evidence that tissue of origin influences gene expression patterns across cancers.

![task_3_hierarchical_clustering_dendrogram](https://github.com/user-attachments/assets/6de39d82-1d12-474f-9b12-2d40f3e3420f)

This dendrogram shows hierarchical clustering of all 9,263 samples based on gene expression profiles. The height of each branch represents the distance (dissimilarity) between clusters. The red horizontal lines indicate the cutoff for 24 clusters (matching the number of cancer types). The complex branching structure reflects the molecular heterogeneity within and between cancer types. While the clusters don't perfectly align with cancer types (7.36% agreement), the major branches likely represent broader biological categories related to tissue origin and cell type.

![task_3_tsne_plot](https://github.com/user-attachments/assets/cc7bcf4e-51c0-4e5c-b2c6-e15b59bad74f)

The t-SNE plot provides a non-linear dimensionality reduction that emphasizes local structure in the data. Each point represents a tumor sample colored by cancer type. The clear separation of samples into distinct clusters demonstrates that gene expression patterns are highly specific to cancer types. Cancer types from similar tissues tend to appear closer in the t-SNE space, suggesting molecular similarities related to tissue of origin. 

![task_3_pca_plot](https://github.com/user-attachments/assets/0958a23d-d2e3-45fb-8421-b355a06a69a1)

Principal Component Analysis (PCA) shows a linear dimensionality reduction of the gene expression data. The first two principal components explain 23.4% of total variance (PC1: 12.6%, PC2: 10.8%). Unlike t-SNE, PCA shows more overlap between cancer types, indicating that the differences between cancer types are complex and not fully captured in two dimensions. Some cancer types form loose clusters, but the boundaries are less distinct than in the t-SNE visualization.
![Task_3_pca_with_kmeans](https://github.com/user-attachments/assets/07d28bb4-ece4-47d0-b1d0-fcf6f2cc4b86)

This plot shows the same PCA projection, but points are colored by k-means cluster assignment (24 clusters) and shaped by cancer type. The limited alignment between colors (clusters) and shapes (cancer types) explains the low agreement score of 5.34%. K-means struggles to separate cancer types in this high-dimensional space, often grouping samples from different cancer types together or splitting single cancer types across multiple clusters. This suggests that standard distance-based clustering may not fully capture the complex biological relationships between cancer types.


#Build a machine learning model to predict cancer type based on gene expression levels.
1. Use Random Forest or Support Vector Machine (SVM)
2. Split data into train/test sets (80/20 split)
3. Evaluate using accuracy, precision, recall, and F1-score

![rf_metrics_by_cancer_type](https://github.com/user-attachments/assets/c08b614e-e5d5-4ba4-8f25-70c9e9054c63)

This chart shows the performance of the Random Forest model for each cancer type across three key metrics: Precision (green), Recall (blue), and F1-Score (red). Most cancer types show excellent performance with metrics above 0.90, indicating the model can accurately identify and distinguish between different cancer types based on gene expression patterns. Cancer types like LAML, DLBC, THCA, and OV show particularly high performance across all metrics (close to 1.0), while cancer types like LUSC and HNSC have slightly lower precision. This suggests that gene expression signatures are highly specific to each cancer type, allowing for accurate classification.

![model_accuracy_comparison](https://github.com/user-attachments/assets/fe8ebb3d-a31e-44b2-b20b-e2a14ff10f19)

This chart compares the overall accuracy of the Random Forest and SVM models. The Random Forest model significantly outperforms SVM with an accuracy of 92.07% versus 31.87% for SVM. This substantial difference suggests that the complex, non-linear patterns in gene expression data are better captured by the ensemble decision tree approach of Random Forest compared to the linear boundaries created by our SVM implementation. This makes Random Forest the preferred model for cancer type prediction from gene expression data.

![rf_confusion_matrix_heatmap](https://github.com/user-attachments/assets/650bdbb5-910b-41a4-910d-92ef60f9de54)

The confusion matrix heatmap visualizes prediction patterns for the top 10 cancer types. The strong diagonal of dark blue squares indicates that most samples are correctly classified into their respective cancer types. The matrix is normalized by row, showing the proportion of each actual cancer type assigned to predicted categories. The near-perfect diagonal pattern with minimal off-diagonal elements demonstrates the model's ability to distinguish between cancer types with high precision. This supports our observation from the t-SNE plot that gene expression patterns are highly specific to cancer types.

These descriptions highlight the key findings from your machine learning analysis:

#Random Forest achieves excellent performance (92% accuracy)

#Different cancer types can be accurately distinguished based on gene expression

#The model performs consistently well across various cancer types

#The confusion matrix confirms minimal misclassification between cancer types

The results complement our earlier clustering analysis by demonstrating that the distinct gene expression patterns observed in the t-SNE plot can be effectively leveraged for automated cancer type classification.


#Use logistic regression or feature importance from Random Forest to find genes that best distinguish one specific cancer type from all others. Which genes have the highest predictive value?

![gene_importance_heatmap](https://github.com/user-attachments/assets/86bf4b63-160f-4c25-ae78-f6749acd3ea2)

This heatmap shows how important each of the top 30 genes is for predicting different cancer types. Darker blue indicates higher importance. Key observations:

TP53 shows high importance across multiple cancer types, particularly for BRCA, LUSC, BLCA, and HNSC. This makes biological sense as TP53 is the most commonly mutated gene in human cancers.
Different cancer types rely on different gene signatures. For example, THCA (thyroid cancer) shows high importance for ESR1 and PPARG, while PRAD (prostate cancer) shows high importance for AR (androgen receptor).
The clustering at the top groups cancer types with similar gene importance patterns. This reveals molecular similarities between cancers (e.g., BRCA and LUAD cluster together).

![top_20_genes_overall](https://github.com/user-attachments/assets/66eb8ed6-b480-4429-a94c-280be1232032)

This bar chart ranks genes by their average importance across all cancer types:

TP53 has the highest overall importance, confirming its central role in cancer biology
ESR1 (estrogen receptor 1) and AR (androgen receptor) rank next, reflecting their roles in hormone-responsive cancers
PPARG is highly ranked, important for metabolism and differentiation
LRRK2, SOX9, and RASSF7 show high importance across multiple cancer types

These genes represent master regulators that influence multiple cancer pathways and can distinguish between different cancer types.

![top_genes_for_ACC](https://github.com/user-attachments/assets/0fe989f2-c4c1-45f5-abb7-9035ce32bb45)

This shows genes specifically important for identifying adrenocortical carcinoma (ACC):

TP53 is again the most important gene
ESR1 and CCND2 (cyclin D2) show high specificity for ACC
Several genes appear important specifically for ACC that aren't in the overall top 20, such as PER2 (a circadian rhythm gene)
The presence of SOX4 and SOX9 (transcription factors) suggests important developmental pathways in ACC.

#Biological Significance
Many of these top genes have established roles in cancer:

#TP53: Tumor suppressor that regulates cell cycle and prevents cancer
#ESR1: Estrogen receptor involved in breast and other hormone-responsive cancers
#AR: Androgen receptor critical in prostate cancer
#CCND2: Cell cycle regulator
#SOX9: Developmental transcription factor increasingly recognized in cancer

The Random Forest model has effectively identified both well-established cancer genes that distinguish between cancer types.









