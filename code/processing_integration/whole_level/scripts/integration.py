import argparse
import json

import anndata as ad
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from anndata import AnnData

sc.logging.print_header()
sc.settings.logging = 4
sc.settings.log_file = 'scanpy.log'
sc.settings.n_jobs = 24

def get_data(data_path: str, meta_path: str, dataset: str):
    adata = sc.read_10x_h5(data_path)
    meta = pd.read_csv(meta_path)
    adata.obs = adata.obs.join(meta)
    adata.obs['dataset'] = dataset
    return adata

def prepare_for_dim_reduction(adata: AnnData):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key = 'sample')
    adata.raw = adata
    var_select = adata.var.highly_variable_nbatches > 5
    var_genes = var_select.index[var_select]
    adata = adata[:, var_genes]
    sc.pp.regress_out(adata, ['nCount_RNA', 'percent.mito'])
    sc.pp.scale(adata, max_value=10)
    return adata

def dim_reduction(adata: AnnData):
    print('DIM REDUCTION STARTS:')
    sc.tl.pca(adata, svd_solver='arpack')
    sce.pp.harmony_integrate(adata, 'sample')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata)
    return adata


def clustering(adata: AnnData, res: list):
    print('CLUSTERING STARTS:')
    for resolution in res:
        sc.tl.leiden(adata, resolution = resolution, key_added = f"leiden_{resolution}")
    return adata


def get_markers(adata: AnnData, res: list):
    for resolution in res:
        sc.tl.rank_genes_groups(adata, groupby = f'leiden_{resolution}', method='wilcoxon',
                                key_added = f'{resolution}', pts = True)
    return adata


def save_adata(adata: AnnData, out_file: str, out_meta: str):
    adata.write_h5ad(out_file)
    adata.obs.to_csv(out_meta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Basic preprocessing: scanpy')
    parser.add_argument('--out_file', type=str, required=True,
                        help='Path to the output h5ad file')
    parser.add_argument('--meta_out_file', type=str, required=True,
                        help='Path to the output adata.obs metadata')
    parser.add_argument('--res', nargs='+', default=[1.0],
                        help='List of the resolutions for the clustering purposes')
    parser.add_argument('--config', type=str, required=True,
                        help='Configuration json file with pathes to h5 / meta data')
    args = parser.parse_args()
    args.res = [float(x) for x in args.res]
    print(args)

    ## READ DATA
    with open(args.config) as file:
        pathes = json.loads(file.read())
    print(pathes)
    adata = ad.concat([get_data(v[0], v[1], k) for k, v in pathes.items()], merge="same")
    adata.obs.to_csv(args.meta_out_file)

    ## NORMALIZATION
    adata = prepare_for_dim_reduction(adata)

    ## DIMENSIONALITY REDUCTION
    adata = dim_reduction(adata)

    ## CLUSTERING
    adata = clustering(adata, args.res)
    adata = get_markers(adata, args.res)

    ## SAVE OUTPUT
    save_adata(adata, args.out_file, args.meta_out_file)

