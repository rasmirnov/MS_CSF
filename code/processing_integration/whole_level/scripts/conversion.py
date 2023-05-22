import argparse
import os
import pandas as pd

import scanpy as sc
from functions.get_description import get_description
from functions.get_expr import get_expr
from functions.get_h5 import get_h5
from functions.get_markers import get_markers
from functions.get_plot_data import get_plot_data


def conversion(adata_path: str, out_dir: str,
                    token: str, name: str, specie: str,
                    link='', description='',
                    public=False, curated=False, max_reduction_dims=5, n_levels=50, use_raw=True):
    species = ['mm', 'hs', 'rn']
    sc.logging.print_header()
    if not specie in species:
        raise ValueError(
            f'Unsupported species - {specie}. \n Supported species are: {", ".join(map(lambda x: x, species))}')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    adata = sc.read_h5ad(adata_path)
    sc.pp.subsample(adata, n_obs=100000)
    if use_raw:
        adata = adata.raw.to_adata()
    adata.obs.columns = adata.obs.columns.str.replace("\\.", "_")
    get_h5(adata=adata, outfile=f'{out_dir}/data.h5')
    get_plot_data(adata, outfile = f'{out_dir}/plot_data.json',
                  max_reduction_dims=max_reduction_dims, n_levels=n_levels)
    get_description(adata=adata, specie=specie,
                    outfile=f'{out_dir}/dataset.json', token=token, name=name,
                    link=link, description=description, public=public, curated=curated)
    get_markers(adata=adata, outfile=f'{out_dir}/markers.json')
    get_expr(adata=adata, outfile=f'{out_dir}/exp_data.json')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert scanpy-processed AnnData to SCN format')
    parser.add_argument('--adata', type=str, required=True,
                        help='Path to h5ad AnnData object preprocessed by scanpy')
    parser.add_argument('--out_dir', type=str, required=True,
                        help='Path to the output directory')
    parser.add_argument('--token', type=str, required=True,
                        help='Unique token for scNavigator')
    parser.add_argument('--name', type=str, required=True,
                        help='Name of the dataset')
    parser.add_argument('--specie', type=str, required=True,
                        help='Species of the dataset. Supported values are: "mm", "hs", and "rn"')
    parser.add_argument('--public', required=False, action='store_true',
                        help='Should be dataset list on the main page in public.')
    parser.add_argument('--description', type=str, required=False, default='', nargs='?',
                        help='Desription of the dataset')
    parser.add_argument('--link', type=str, required=False, default='', nargs='?',
                        help='External link for the dataset')
    parser.add_argument('--curated', required=False, action='store_true',
                        help='Should be dataset list on the main page in curated')
    parser.add_argument('--max_reduction_dims', type=int, required=False, default=5, nargs='?',
                        help='Should be dataset list on the main page in curated')
    parser.add_argument('--n_levels', type=int, required=False, default=300, nargs='?',
                        help='Specify in case you have more levels in your target factor variable than 50'
                             '(e.g., more than 50 clusters in the certain resolution)')
    args = parser.parse_args()
    print(args)
    conversion(adata_path=args.adata, out_dir=args.out_dir, token=args.token, name=args.name, specie=args.specie,
               link=args.link, description=args.description, public=args.public,curated=args.curated,
               max_reduction_dims=args.max_reduction_dims,n_levels=args.n_levels)
