import numpy as np
import os
import anndata as ad
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from anndata import AnnData
import matplotlib as plt

microglia = sc.read_h5ad(filename = 'path to your subsetted object')
adata = microglia
adata = adata.raw.to_adata()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, 
min_disp=0.5, batch_key = 'dataset')
adata.raw = adata
var_select = adata.var.highly_variable_nbatches > 2
var_genes = var_select.index[var_select]
adata = adata[:, var_genes]
sc.pp.regress_out(adata, ['nCount_RNA', 'percent.mito'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sce.pp.harmony_integrate(adata, 'study')  
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=5, n_pcs=50)
sc.tl.umap(adata)
res = [0.1, 0.2, 0.3,
       0.4, 0.5, 0.6]

for resolution in res:
    sc.tl.leiden(adata, resolution = resolution, key_added = f"leiden_{resolution}")
# gene markers
for resolution in res:
    sc.tl.rank_genes_groups(adata, 
                            groupby = f'leiden_{resolution}',
                            method='wilcoxon',
                            key_added = f'{resolution}',
                            pts = True)
# save
adata.write_h5ad(filename =  'destination path with filename')