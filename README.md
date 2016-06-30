# GiniClust

GiniClust is a clustering method implemented in R for detecting rare cell types from large scale single cell gene expression data. It can be applied to data set come from different platform, such as multiplex qPCR data, tranditional single cell RNAseq or newly emerging UMI-based single cell RNAseq, such as inDrops and Drop-seq. GiniClust is created and maintained by the Yuan Lab at DFCI.


# Preparation
Download the zip file GiniClust.zip and unzip it.  The unzipped folder will contain one script 'Giniclust_Main.R', one folder 'Rfunction' including all functional R codes and two data samples 'Data_GBM.csv' and 'Data_qPCR.csv'.

Make sure 'Giniclust_Main.R', 'Rfunction' and input data in the same directory.

# Input File Format
The input file is a gene expression matrix in csv file format.
Especially, for qPCR data, each entry is log2 gene expression level and for RNA-seq data, each entry is UMI-Count/Cell or Raw-Read-Count/Cell.  The first row contains cell names. The first column contains unique gene ids. 

For example, in R 
```R
>ExprM.RawCounts  <- read.csv("Data_GBM.csv", sep=",", head=T)
>ExprM.RawCounts[1:4,1:4]
```
you can take a look at one of our test dataset:

|Table   |MGH26 | MGH26.1 | MGH26.2 | MGH26.3|
|------------ |------------| -------------|------------ | -------------|
|1/2-SBSRNA4| 0      |47       |0       |0|
|A1BG          | 41      |80       |3       |0|
|A1BG-AS1        |0       |0       |0      |0|
|A1CF            |0       |0       |0       |0|


# Rfunction
- GiniClust_parameters.R: All the parameters used in Giniclust
- GiniClust_packages.R: Installing and loading required packages
- GiniClust_Preprocess.R: Preprocessing the input data
- GiniClust_Filtering.R: Filtering preprocessed data
- GiniClust_Fitting.R: LOESS curve fitting and selecting final genes which are used for clustering
- GiniClust_Clustering.R: DBSCAN for clustering
- GiniClust_tSNE.R: Visulization using tsne
- DE_MAST.R: MAST analysis for RNA seq, which is used for identifing differential genes between rare cell types and other cells
- DE_MAST_figures.R:  generating MAST figures
- DE_t_test.R: t-test for qPCR data


# Usage 
    Rscript Giniclust_Main.R [options]

Options:
- -f CHARACTER or  --file=CHARACTER, input dataset file name 
- -t CHARACTER or --type=CHARACTER, input dataset type: choose from 'qPCR' or 'RNA-seq' 
- -o CHARACTER or  --out=CHARACTER, output folder name [default= results]
- -h or  --help, Show help message and exit

For example, running GiniClust with 'Data_GBM.csv',
```sh
$ Rscript GiniClust_Main.R -f Data_GBM.csv -t RNA-seq -o GBM_results
```
running GiniClust with 'Data_qPCR.csv',
```sh
$ Rscript GiniClust_Main.R -f Data_qPCR.csv -t qPCR -o  qPCR_results
```


# Results
Two folders will be created automatically in the currently directory. 
One is the output folder with a default or specified name, which includes all figures and tables in our manuscript.
output file list and description:

- Dataname_rawCounts.csv: the raw counts 
- Dataname_normCounts.csv: the normalized counts
- Dataname_gene.expression.matrix.RawCounts.filtered.csv: the raw counts after filtering 
- Dataname_gene.expression.matrix.normCounts.filtered.csv: the normalized counts after filtering  
- Gini_related_table_RNA-seq.csv: the table related with Gini index for RNA-seq data
- Gini_related_table_qPCR.csv: the table related with Gini index for qPCR data
- Dataname_clusterID.csv: clustering result, the first column represents cell IDs and the second column is the corresponding cluster result for each cell.
- Dataname_Rtnse_coord2.csv: coordinates of cells in tSNE plot 
- Dataname_bi-directional.GiniIndexTable.csv: For qPCR data the table of bidirectional Gini index
- RareCluster_lrTest.csv: lrTest results by MAST analysis for RNA-seq data 
- RareCluster.diff.gene.t-test.results.csv t-test results for qPCR data
sub-folder 'figures':
- Dataname_histogram of Normalized.Gini.Socre.pdf: histogram of estimated p-values based on a normal distribution approximation for genes
- Dataname_smoothScatter_pvalue_gene.pdf: the smoothScatter plot in which the red points are the selected high Gini genes according to specified cutoff
- Dataname_tsne_plot.pdf: tSNE plot for cells
- Dataname_RareCluster_diff_gene_overlap.pdf: Venn diagram for differentially expressed genes and high gini genes
- Dataname_RareCluster_overlapgene_rawCounts_bar_plot.genename.pdf: barplot of rare cluster and major cluster for the overlap genes

The other folder is the folder 'Libarary' , which includes all newly installed packages


# Publication 
Jiang L, Chen H, Pinello L, Yuan GC. GiniClust: Detecting rare cell types from single-cell gene expression data with Gini Index. Genome Biology (in press)

# Contact 
Lan Jiang ( lan_jiang at hms dot harvard dot edu ) or Guo-Cheng Yuan ( gcyuan at jimmy dot harvard dot edu )

# License
The source code is released under the MIT license (https://opensource.org/licenses/MIT).

# Troubleshooting
For MAC and WINDOWS only the official R installation file is supported and tested. Using other installation methods, such as brew, may lead to running error.

Some users have encountered errors due to problems with MAST installation. Please refer to the MAST package URL (https://github.com/RGLab/MAST) for help. 
