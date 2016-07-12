The present README.txt file provides a brief overview of the R scripts behind much of GiniClust's computations.

- GiniClust_parameters.R: Parameter setting for Giniclust
- GiniClust_packages.R: Installation, includiong loading additional R packages required for GiniClust
- GiniClust_Preprocess.R: Preprocessing the input data
- GiniClust_Filtering.R: Filtering the preprocessed data
- GiniClust_Fitting.R: Normalizing the Gini index by LOESS curve fitting, selecting a subset of the genes for clustering
- GiniClust_Clustering.R: Clustering the cells based on the selected genes using DBSCAN
- GiniClust_tSNE.R: Visualization of the clustering results using t-SNE
- DE_MAST.R: Differential gene expression analysis for RNAseq data using MAST
- DE_MAST_figures.R: Generating MAST-like figures
- DE_t_test.R: Differential gene expression analysis for qPCR data using t-test
