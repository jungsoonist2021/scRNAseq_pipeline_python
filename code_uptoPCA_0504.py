# -*- coding: utf-8 -*-
"""

by Angelo Digirolamo & Jung Soo jason Kim

"""

#C:/Users/jungs/10x/matrix.mtx

import sklearn 
import matplotlib.pyplot as plt
import pandas as pd
import csv
import gzip
import os
import scipy.io

from scipy.io import mmread
import numpy as np


from gtfparse import read_gtf


#this is to get all specific neuronal subtypes

matrix_dir = "C:/Users/jungs/c_elegans"

features_path_2 = os.path.join(matrix_dir, "cell_type_annotation_lookup_table.csv.gz")

#read csv into dataframe
cell_type = [row[0] for row in csv.reader(gzip.open(features_path_2, mode="rt"), delimiter="\t")]

#%%

#convert list into dataframe
cell_type2 = pd.DataFrame(cell_type)

#%%

#convert first row (index 0) into header/column
cell_type2.columns = cell_type2.iloc[0]
cell_type2 = cell_type2[1:]

#%%

#add new column 'list' that splits string by comma

cell_type2['list'] = cell_type2['barcode,cell.type,tissue.type'].str.split(',')

#%%

#split 'list'column into separate columns
cell_type3 = pd.DataFrame(cell_type2['list'].tolist()).fillna('')

#%%

#get rows that have 'Neuron' in Column 2
cell_type4 = cell_type3.loc[cell_type3[2] == 'Neuron']



neurons = cell_type4[0]

#%%
#acquire all neuronal subtypes: 130 

cell_type_neurons = set(cell_type4[1])

neuron_type_list = list(cell_type_neurons)


#%%
#acuire raw anndata set from Taylor et al
import scipy.io

from scipy.io import mmread

import scanpy as sc
import os


dir1 = "C:/Users/jungs/c_elegans/taylor2020.h5ad"


adata2 = sc.read(dir1)


#%%
#filter anndata to focus on neuronal data
adata_neurons= adata2[adata2.obs.cell_type.isin(neuron_type_list)]

#%%
#filter out genes expressed in less than 5 cells

sc.pp.filter_genes(adata_neurons, min_cells=5)




#%%
# preprocessing (log normalization, scaling) before computing PC

#remove non-neuronal genes??

adata_cpm = adata_neurons.copy() # apply this to a copy so we can compare methods
adata_cpm.raw = adata_cpm # store a copy of the raw values before normalizing
sc.pp.normalize_per_cell(adata_cpm, counts_per_cell_after=1e6)
sc.pp.log1p(adata_cpm)
sc.pp.scale(adata_cpm)

#%%
#Batch Correction!
#sc.external.pp.mnn_correct(adata_cpm,batch_key='sample_batch',k=5)



#%%
#run PCA
adata_pca2 = sc.pp.pca(adata_cpm,n_comps = 125)
#sc.pp.neighbors(adata_cpm, n_neighbors=75)

#125 PCs, umap.min_dist =0.3, umap.n_neighbors = 75, alignment_k (for align_cds) = 5, clustering resolution 3e-3).


#%%

adata_cpm.write('anndata_cpm_afterPCA.h5ad')





