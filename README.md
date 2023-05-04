# scRNAseq_pipeline_python
Recreating results from "Molecular Topography of an Entire Nervous System" (Taylor et al., 2021) using Python pipeline

The paper can be found here : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136049 
The associated R program used can be found here : https://github.com/cengenproject/Initial_single_cell_analysis

The purpose of the program is to label cell identities represented within a scRNA seq dataset. Furthermore, this program is designed to use those generated cell identities to generate dot plots that help visualize expression patterns of genes of interest. This program assembles a range of pre-exsiting python functions with the intent to match some capabilities of the R program Monocle 3. This program is flexible enough to adapt to other scRNA seq datasets with some modification of the code.

The data processing pipeline is broken down as follows : 

1. Obtaining the dataset as an anndata object 

2. Preprocessing: Filtering out lowly expressed genes 

3. Batch Correction 

4. PCA 

5. Dimensionality Reduction with UMAP

6. Leiden algorithm for community detection 

7. Sub-clustering of Centroid partition: PAGA 

8. scoreCT 

9. Generation of dot plots for gene subsets of interest 

10. Cluster Purity analysis 


Hardware Requirements:
This program was run using a Samsung DESKTOP 2PCJ0FG with 16 GB of RAM. It is not known whether devices with less RAM can successfully run the program. 

This program uses the following python libraries: 
1. Sklearn
2. Matplotlib
3. Pandas
4. Csv
5. Gzip 
6. Os
7. Scipy.io 
8. Numpy
9. Seaborn
10. Scanpy
11. scoreCT






