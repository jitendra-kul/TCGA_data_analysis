# TCGA_data_analysis

Questions:
Load the TCGA dataset and provide a summary including:
   1. The number of samples per cancer type
   2. The mean and standard deviation of expression levels for each gene
   3. The top 5 most variable genes across all samples

Another question is"Identify the top 10 differentially expressed genes between 5 cancer types using a Kruskal-Wallis test test. Visualize the expression of these genes with a heatmap. (you can make multiple visualization plots)"
![gene_correlation_heatmap](https://github.com/user-attachments/assets/5462f889-4215-4505-9511-a00b7fb325b3)
![avg_expression_heatmap](https://github.com/user-attachments/assets/26965a4c-e45c-4bee-8f4b-d71c866043a8)
![top_genes_heatmap](https://github.com/user-attachments/assets/83da5f0a-6b86-4977-a0b6-cd0d193299ee)

Perform hierarchical clustering or k-means clustering on the dataset using all genes.
   
•	Visualize the clustering with a dendrogram or t-SNE/PCA plot.

•	Interpret whether cancers with similar tissue origins cluster together


Here I performed clustering to determine if cancers with similar tissue origins have similar gene expression patterns.

Overview of the dataset: From the summary file which will be generated after a user will run the task3 code: it shows 9263 samples, 716 genes and 24 cancer types.

![Task_3_cancer_correlation_heatmap](https://github.com/user-attachments/assets/a1b156b7-69c9-446e-a026-f1c8d0b8760a)
This heatmap shows correlations between different cancer types based on their gene expression profiles. Red indicates high correlation (similarity) while blue indicates low correlation. The hierarchical clustering (dendrogram) at the top groups cancer types with similar expression patterns. Note how cancers from similar tissue origins cluster together: kidney cancers (KIRC, KIRP, KICH), brain tumors (GBM, LGG), and squamous cell carcinomas (BLCA, LUSC, HNSC, CESC) form distinct groups. This provides strong evidence that tissue of origin influences gene expression patterns across cancers.

![task_3_hierarchical_clustering_dendrogram](https://github.com/user-attachments/assets/ac18a452-076b-402f-9ca6-e3235b77a033)
This dendrogram shows hierarchical clustering of all 9,263 samples based on gene expression profiles. The height of each branch represents the distance (dissimilarity) between clusters. The red horizontal lines indicate the cutoff for 24 clusters (matching the number of cancer types). The complex branching structure reflects the molecular heterogeneity within and between cancer types. While the clusters don't perfectly align with cancer types (7.36% agreement), the major branches likely represent broader biological categories related to tissue origin and cell type.

![task_3_tsne_plot](https://github.com/user-attachments/assets/3403ba46-ee4a-4e09-ae40-f28f05c6bf43)
The t-SNE plot provides a non-linear dimensionality reduction that emphasizes local structure in the data. Each point represents a tumor sample colored by cancer type. The clear separation of samples into distinct clusters demonstrates that gene expression patterns are highly specific to cancer types. Cancer types from similar tissues tend to appear closer in the t-SNE space, suggesting molecular similarities related to tissue of origin. 

![task_3_pca_plot](https://github.com/user-attachments/assets/3b78607d-40b4-4925-8e0e-a3ccc5389d22)
Principal Component Analysis (PCA) shows a linear dimensionality reduction of the gene expression data. The first two principal components explain 23.4% of total variance (PC1: 12.6%, PC2: 10.8%). Unlike t-SNE, PCA shows more overlap between cancer types, indicating that the differences between cancer types are complex and not fully captured in two dimensions. Some cancer types form loose clusters, but the boundaries are less distinct than in the t-SNE visualization.

![Task_3_pca_with_kmeans](https://github.com/user-attachments/assets/478c9fe1-50a8-4f02-a9e2-542e8c605552)
This plot shows the same PCA projection, but points are colored by k-means cluster assignment (24 clusters) and shaped by cancer type. The limited alignment between colors (clusters) and shapes (cancer types) explains the low agreement score of 5.34%. K-means struggles to separate cancer types in this high-dimensional space, often grouping samples from different cancer types together or splitting single cancer types across multiple clusters. This suggests that standard distance-based clustering may not fully capture the complex biological relationships between cancer types.


