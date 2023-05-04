# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 22:38:36 2023

@author: jungs
"""

import scipy.io

from scipy.io import mmread

import scanpy as sc
import os


dir5 = "C:/Users/jungs/c_elegans/anndata_cpm_afterPCA2.h5ad"


adatawpca = sc.read(dir5)

#%%
# looked at pca components


sc.pl.pca(adatawpca)
sc.pl.pca_overview(adatawpca)


#%%
#run umap preprocessing with default parameters

sc.external.pp.mnn_correct(adatawpca,batch_key='sample_batch',k=5)

sc.pp.neighbors(adatawpca, n_neighbors=75)
sc.tl.umap(adatawpca, min_dist=0.3, spread=1, random_state=1, n_components= 2)


#%%
#run leiden clustering algoirthm
sc.tl.leiden(adatawpca, resolution = 1.5)
sc.pl.umap(adatawpca, color='leiden')
sc.pl.umap(adatawpca, color='leiden', legend_loc = 'on data')


#we get 88 initial clusters using resolution 1.5



#%%
sc.pl.pca(adatawpca, color="leiden")

sc.pl.pca_variance_ratio(adatawpca)
sc.pl.pca_loadings(adatawpca, include_lowest=False)


#%%
#check purity of initial clusters

import pandas as pd
#mean number of cells per cluster

tmp = pd.crosstab(adatawpca.obs['leiden'],adatawpca.obs['cell_type'], normalize='index')
tmp.plot.bar(stacked=True).legend(loc='upper right')



#%%
tmp.to_csv('tmp.csv',index=True)

#%%
#we settled at resolution 1.5 because we were left with three clusters with 100% purity 
initial_pure_count = tmp.T[tmp.T==1].count()

#to verify cell-type composition of cluster 84/85/87
cluster87 = adatawpca[adatawpca.obs['leiden'].isin(['87']),:]

print(cluster87.obs['cell_type'])



#%%
# ran paga to focus on centroid group of clusters

sc.tl.paga(adatawpca, groups = 'leiden')

sc.pl.paga(adatawpca, threshold=1, color='leiden', fontsize=10, edge_width_scale=0.1)




#%%
#focused on 27 clusters with high connectivity
#adata_sub2 = adatawpca[adatawpca.obs['leiden'].isin(['29','32','83','64','62','4','3','86','45','70','44','87','26','13','53','85','33','46','52','51','82','58','73','0','84'])]

adata_sub2 = adatawpca[adatawpca.obs['leiden'].isin(['67','30','64','41','65','85','0','51','87','86','54','81','82','79','45','55','18','42','17','6','2','80','14','4','62','39','37'])]

#%%
#preprocessing
sc.pp.filter_genes(adata_sub2, min_cells=0)


adata_pca3 = sc.pp.pca(adata_sub2, n_comps = 110)

sc.external.pp.mnn_correct(adata_sub2,batch_key='sample_batch',k=50)

sc.pp.neighbors(adata_sub2,n_neighbors = 50, n_pcs=110)
#sc.external.pp.mnn_correct(adata_cpm,k=50)


sc.tl.umap(adata_sub2, min_dist=0.3, spread=1.0, random_state=1, n_components= 2)

#%%
# leiden clustering and umap

sc.tl.leiden(adata_sub2, resolution = 1)
sc.pl.umap(adata_sub2, color='leiden') #32 clusters from 11 cluster subset (old)
                                        #44 clusters from 20 clusters (1.25)
                                        #45 clusters from 25 clusters (min_dist = 0.2!!!)
                                        #48 clusters from 27 clusters???

#88 -27+ 48 = 109 clusters
#%%
#added "0_" to clusters that were chosen to be sub-clustered
adata_sub2.obs["leiden"] = "0_" + adata_sub2.obs["leiden"].astype(str)
#%%
#visualization check
sc.pl.umap(adata_sub2, color='leiden')

#%%
#newly subclustered clusters are updated into main list of leiden clusters
adatawpca.obs["leiden"] = adatawpca.obs["leiden"].astype(str)
adatawpca.obs["leiden"].update(adata_sub2.obs["leiden"])

sc.pl.umap(adatawpca, color='leiden')

#%%
#find cell_type composition of all 110 clusters
import pandas as pd


tmp2 = pd.crosstab(adatawpca.obs['leiden'],adatawpca.obs['cell_type'], normalize='index')
tmp2.plot.bar(stacked=True).legend(loc='upper right')

tmp2.to_csv('tmp2.csv',index=True)
    

##find how many clusters are pure, how many satisfy the 99% threhold


#%%
##find how many clusters are pure, how many satisfy the 99% threhold


result = tmp2.dtypes
print(result)

tmp2_new = tmp2.astype(int)

tmp2_transposed=tmp2_new.T

#c = tmp2.T.value_counts()['1']
count1 = tmp2_transposed[tmp2_transposed==1].count()
print(count1.value_counts()) #this shows there are 'five' 100% pure clusters out of 110


#%%

tmp2_float_transposed = tmp2.T

count2 = tmp2_float_transposed[tmp2_float_transposed>0.9900000000].count()

print(count2.value_counts()) #this shows there are 26 clusters with >99% purity


#%%
#this shows the number of cells per each cluster
cluster_counts = adatawpca.obs['leiden'].value_counts()


#%%
#get average number of cells per cluster
avg = cluster_counts.mean()


#%%
#create anndata file with full leiden cluster annotation

adatawpca.write('anndata_withsubclusters_1.5.h5ad')


#%%


dir = "C:/Users/jungs/c_elegans/anndata_withsubclusters_1.5.h5ad"


adatawpca = sc.read(dir)




#%%
#install necessary packages to run scoreCT
import os
import sys
import pandas as pd
import scanpy as sc
import scorect as ct
# Ignore warnings in this tutorial
import warnings
warnings.filterwarnings('ignore')

#%%
#run scanpy rank genes function for all leiden clusters
adatawpca.uns['log1p']["base"] = None

sc.tl.rank_genes_groups(adatawpca, groupby='leiden', n_genes=len(adatawpca.raw.var), use_raw=True)

#%%
# Use scoreCT one line API for assignments
ct.scorect(adatawpca,
          marker_path="C:/Users/jungs/c_elegans/WBG2.gmt",
          K_top=300,
          m_bins=5,
          null_model='multinomial',
          cluster_key='leiden')


#%%
adatawpca.write('anndata_usingShortList_1.5_scoreCT.h5ad')

#adata_withsubclusters.write('anndata_afterscoreCT.h5ad')


#%%
dir = "C:/Users/jungs/c_elegans/anndata_usingShortList_1.5_scoreCT.h5ad"


adatawpca = sc.read(dir)

#%%
sc.pl.umap(adatawpca, color='leiden')

#%%
#plot umap of correct cell_type annotation per cluster compared to scoreCT predicted result 
sc.pl.umap(adatawpca, color=['cell_type', 'scoreCT'], title=['True','Predicted'], wspace=0.33, legend_loc = 'on data')


#sc.pl.umap(adatawpca, color=['scoreCT'], title=['Predicted'], wspace=0.33, legend_loc = 'on data')


#%%
#focus on cells/clusters that scoreCT successfully labeled
subset_match_neighbor = adatawpca[adatawpca.obs['scoreCT'].astype(str)==adatawpca.obs['cell_type'].astype(str)]



#%%
#get number of correctly labeled cell_types
print(subset_match_neighbor.obs['scoreCT']) # 63 /55 (resolutin 1.5)/ 44

#%%
#number of cell_types that scoreCT identified
print(adatawpca.obs['scoreCT']) # 77


#%%
#number of cell_types in original anndataset including 10 subtypes 
print(adatawpca.obs['cell_type']) # 130/118



#%%
#number of cells per correctly labeled cluster
cluster_counts2 = subset_match_neighbor.obs['scoreCT'].value_counts()

#average number of cells for correctly labeled clusters
avg2 = cluster_counts2.mean()

#%%

#69 clusters, 63 cell types: there are some clusters that were labeled the same..
#clusters (16,38),(AS1)
# (0_13, 15), (RMF)
#(0_1, 0_11) (SIA)
# (0_5, 0_31) (SIB)
#(0_20,0_25) (URA)
#(1,14) (VC)

#cell_type composition analysis for correctly labeled clusters

tmp3 = pd.crosstab(subset_match_neighbor.obs['leiden'],subset_match_neighbor.obs['cell_type'], normalize='index')
tmp3.plot.bar(stacked=True).legend(loc='upper right')

tmp3.to_csv('tmp3.csv',index=True)
    
#%%
#get cluster names: tmp3 from correctly labeled clusters, tmp2 is from all leiden clusters
tmp3_clusters= list(tmp3.index)

tmp2_clusters = list(tmp2.index)


#%%
#just look at leiden clusters that were correctly labeled, and count number of 100% pure clusters
tmp2_filtered = tmp2.filter(items = tmp3_clusters, axis=0)

tmp2_filtered_transposed = tmp2_filtered.T

count3 = tmp2_filtered_transposed[tmp2_filtered_transposed==1].count()

print(count3.value_counts()) # 5 100% pure clusters

#%%

 #this gives number of clusters with 99% purity that were correctly labeled
count4 = tmp2_filtered_transposed[tmp2_filtered_transposed>0.9900000].count()
             

print(count4.value_counts()) #17 clusters that were above 99percent purity     

             

#%%

homeobox_gene_name_list = ["ceh-13","lin-39","mab-5","egl-5","nob-1","php-3","pal-1","vab-7","ceh-1","ceh-9","ceh-19","ceh-22","ceh-24","ceh-27","ceh-28","ceh-30","ceh-31","ceh-51","cog-1","mls-2","tab-1","vab-15","ceh-2","ceh-5","ceh-7","ceh-12","ceh-16","ceh-23","ceh-43","ceh-62","ceh-63","pha-2","ceh-21","ceh-39","ceh-41","ceh-38","ceh-44","ceh-48","ceh-49","dve-1","ceh-14","lin-11","lim-4","lim-6","lim-7","mec-3","ttx-3","ceh-6","ceh-18","unc-86","eyg-1","pax-3","vab-3","alr-1","ceh-8","ceh-10","ceh-17","ceh-36","ceh-37","ceh-45","ceh-53","ceh-54","dsc-1","unc-30","unc-4","unc-42","ttx-1","pros-1","hmbx-1","ceh-32","ceh-33","ceh-34","unc-39","ceh-20","ceh-40","ceh-60","irx-1","unc-62","zag-1","zfh-2","ceh-58","ceh-75","ceh-79","ceh-87","ceh-88","ceh-92","ceh-93","nsy-7","ceh-57","ceh-74","ceh-76","ceh-81","ceh-82","ceh-83","ceh-84","ceh-85","ceh-86","ceh-89","ceh-90","ceh-91","ceh-99","ceh-100","duxl-1"]
homeobox_wormbase_identifier_list_full = ["WBGene00000437","WBGene00003024","WBGene00003102","WBGene00001174","WBGene00003779","WBGene00004024","WBGene00003912","WBGene00006873","WBGene00000428","WBGene00000434","WBGene00000442","WBGene00000445","WBGene00000447","WBGene00000449","WBGene00000450","WBGene00000451","WBGene00000452","WBGene00013583","WBGene00000584","WBGene00003377","WBGene00006380","WBGene00006881","WBGene00000429","WBGene00000430","WBGene00000432","WBGene00000436","WBGene00000439","WBGene00000446","WBGene00000463","WBGene00011069","WBGene00045215","WBGene00004011","WBGene00000444","WBGene00000460","WBGene00000462","WBGene00000459","WBGene00000464","WBGene00015934","WBGene00017538","WBGene00022861","WBGene00000438","WBGene00003000","WBGene00002987","WBGene00002988","WBGene00002989","WBGene00003167","WBGene00006654","WBGene00000431","WBGene00000441","WBGene00006818","WBGene00013147","WBGene00003939","WBGene00006870","WBGene00044330","WBGene00000433","WBGene00000435","WBGene00000440","WBGene00000457","WBGene00000458","WBGene00022837","WBGene00015651","WBGene00020485","WBGene00001096","WBGene00006766","WBGene00006744","WBGene00006778","WBGene00006652","WBGene00000448","WBGene00018786","WBGene00000453","WBGene00000454","WBGene00000455","WBGene00006775","WBGene00000443","WBGene00000461","WBGene00017690","WBGene00007984","WBGene00006796","WBGene00006970","WBGene00022518","WBGene00007417","WBGene00008242","WBGene00007749","WBGene00018022","WBGene00008195","WBGene00013431","WBGene00019864","WBGene00044508","WBGene00007416","WBGene00013876","WBGene00022395","WBGene00018434","WBGene00018433","WBGene00018446","WBGene00016557","WBGene00019137","WBGene00018355","WBGene00009231","WBGene00010995","WBGene00013425","WBGene00044032","WBGene00012584","WBGene00022554"]

homeobox_genes_not_found = ['WBGene00000432', 'WBGene00000461', 'WBGene00008242', 'WBGene00013147', 'WBGene00018433', 'WBGene00019137', 'WBGene00022395']

homeobox_wormbase_identifier_list = list(set(homeobox_wormbase_identifier_list_full) - set(homeobox_genes_not_found))

homeobox_gene_dict = {}

'''
for i in range(len(homeobox_gene_name_list)):
    for j in range(len(homeobox_wormbase_identifier_list)):
        homeobox_gene_dict[homeobox_wormbase_identifier_list[i]]= homeobox_gene_name_list[j]

print(homeobox_gene_dict)
'''


#%%
neuropeptide_wormbase_identifier_list_full = ["WBGene00000085","WBGene00016090","WBGene00016748","WBGene00015005","WBGene00007146","WBGene00015299","WBGene00015416","WBGene00007483","WBGene00015921","WBGene00016044","WBGene00008142","WBGene00016858","WBGene00008300","WBGene00020712","WBGene00021439","WBGene00017015","WBGene00010191","WBGene00022605","WBGene00019262","WBGene00020524","WBGene00020523","WBGene00020525","WBGene00194821","WBGene00021275","WBGene00012993","WBGene00017038","WBGene00012962","WBGene00013222","WBGene00016428","WBGene00007951","WBGene00022606","WBGene00008481","WBGene00001175","WBGene00017433","WBGene00008885","WBGene00009333","WBGene00018089","WBGene00018222","WBGene00009931","WBGene00010058","WBGene00018924","WBGene00010143","WBGene00010179","WBGene00010315","WBGene00010329","WBGene00015323","WBGene00019019","WBGene00019444","WBGene00019445","WBGene00019448","WBGene00019496","WBGene00010735","WBGene00020023","WBGene00011765","WBGene00020586","WBGene00021510","WBGene00007346","WBGene00016149","WBGene00016909","WBGene00008342","WBGene00017661","WBGene00018191","WBGene00018728","WBGene00009965","WBGene00018798","WBGene00007614","WBGene00013871","WBGene00016570","WBGene00019224","WBGene00008736","WBGene00008737","WBGene00013642","WBGene00010375","WBGene00019231","WBGene00019367","WBGene00019370","WBGene00002252","WBGene00019774","WBGene00016747","WBGene00019616","WBGene00017176","WBGene00003807","WBGene00008278","WBGene00016110","WBGene00020689","WBGene00013883","WBGene00012275","WBGene00012084","WBGene00006864","WBGene00015559","WBGene00008065","WBGene00015364","WBGene00003808","WBGene00011578","WBGene00020727","WBGene00022004","WBGene00021328","WBGene00020086","WBGene00011381","WBGene00011372","WBGene00018344","WBGene00018886","WBGene00013848","WBGene00007006","WBGene00019183","WBGene00020319","WBGene00013782","WBGene00009278","WBGene00013187","WBGene00016842","WBGene00007635","WBGene00021983","WBGene00009619","WBGene00018067","WBGene00016984","WBGene00013974","WBGene00011582","WBGene00008808","WBGene00015735","WBGene00020008","WBGene00020069","WBGene00014035","WBGene00007664","WBGene00010986","WBGene00009629","WBGene00022086","WBGene00005094","WBGene00015401","WBGene00013375","WBGene00005862","WBGene00005865","WBGene00005869","WBGene00005870","WBGene00020138","WBGene00011371","WBGene00020420","WBGene00011709","WBGene00020664","WBGene00006428","WBGene00016265","WBGene00021124","WBGene00021357","WBGene00021497","WBGene00021524","WBGene00013496","WBGene00014248","WBGene00022824","WBGene00014122","WBGene00003862","WBGene00006576","WBGene00003863","WBGene00016761"]

neuropeptide_genes_not_found = ["WBGene00003862", "WBGene00003863", "WBGene00005862", "WBGene00005865", "WBGene00005869", "WBGene00005870", "WBGene00020525"]

neuropeptide_wormbase_identifier_list = list(set(neuropeptide_wormbase_identifier_list_full) - set(neuropeptide_genes_not_found))

sc.pl.dotplot(
    subset0,
    groupby="scoreCT",
    var_names = homeobox_wormbase_identifier_list,
    standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
    title = "Mean Expression of Homeobox Genes per Cluster",
    swap_axes = True,
    save = "Mean Expression of Homeobox Genes per Cluster Dotplot.png"
)

sc.pl.matrixplot(subset0, var_names = homeobox_wormbase_identifier_list, groupby ='scoreCT', use_raw=False, vmax = 6,title = "Mean Expression of Homeobox Genes per Cluster",swap_axes = True,save = "Mean Expression of Homeobox Genes per Cluster Matrixplot.png")







