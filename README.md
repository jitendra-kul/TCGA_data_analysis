# TCGA_data_analysis

Questions:
Load the TCGA dataset and provide a summary including:
   1. The number of samples per cancer type
   2. The mean and standard deviation of expression levels for each gene
   3. The top 5 most variable genes across all samples

Another question is"Identify the top 10 differentially expressed genes between 5 cancer types using a Kruskal-Wallis test test. Visualize the expression of these genes with a heatmap. (you can make multiple visualization plots)"

Another question is"Identify the top 10 differentially expressed genes between 5 cancer types using a Kruskal-Wallis test test. Visualize the expression of these genes with a heatmap. (you can make multiple visualization plots)"
![gene_correlation_heatmap](https://github.com/user-attachments/assets/478480a0-68b0-411c-a198-c67fc68c0af4)
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








