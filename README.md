# enrichment_analysis_celltype
Cell type enrichment analysis using gene signatures and cluster markers

# Functions

## enrich

The **enrich** function performs cell type enrichment analysis using signatures from annotation databases and cluster markers. It first computes a contingency table for every pair of comparison between cell-type annotations and cluster IDs and then performs a Fisher's exact test of independence.

### Input parameters

The **enrich** function takes as input 5 parameters:

- clusters - path of the file from which the cluster information are to be read in. The cluster dataset should be a csv file with genes along rows and columns corresponding to adjusted p value, log fold change (or average difference) and cluster ID respectively. An example file can be loaded from the data folder.
- pval - adjusted p value threshold to select significant gene memberships for clusters
- lfc - log fold change threshold to select significant gene memberships for clusters
- annoDB - annotation database to be used for cell-type signatures. One of 'xCell' or 'brain_anatomy' can be selected. 'xCell' enumerates 64 immune and stromal cell-types using a novel gene-signature based method. 'brain_anatomy', on the other hand, applies a correlation-based metric called differential stability to identify 21 significant gene expression signatures representing distinct anatomical regions of the brain (default = 'xCell').
- plot - a logical value indicating whether to plot the top 3 cluster cell-types as a barplot (default = TRUE)

### Usage of the enrich function

The **enrich** function will display a barplot (if plot = TRUE) of the top 3 significantly enriched cluster cell-types from the Fisher's exact test and also returns a list containing the following components:
- annotation_adj_matrix - adjacency matrix containing genes (rows) and cell-type annotations (columns) used in the Fisher's test
- clusters - clusters dataset with genes along rows and adjusted p value, average log fold change (or average difference) and cluster ID along columns
- matrix_fisher - the results of the Fisher's test containing the p-value of the test, an estimate of the odds ratio, the number of genes present in the corresponding signature, the number of genes present in the cluster, genes common to both and the corresponding cluster ID. 
- significant_fisher - top 3 significantly enriched results from the Fisher's test also used for making the barplot

### Example
