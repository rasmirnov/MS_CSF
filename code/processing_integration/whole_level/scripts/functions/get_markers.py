import json
import re

import numpy as np
import scanpy as sc
from anndata import AnnData


def get_markers(adata: AnnData, outfile: str):
    resolutions = [x for x in adata.obs.select_dtypes(include='category').columns if re.findall('leiden|louvain', x)]
    markers = dict()
    for res in resolutions:
        res_key = re.sub('_', '.', re.findall('[0-9].*[0-9]', res)[0])
        markers_ta = sc.get.rank_genes_groups_df(adata, key=res_key,
                                                 group=list(adata.uns[res_key]['names'].dtype.names),
                                                 pval_cutoff=0.01).round(4)
        markers_ta.rename(columns={'group': 'cluster', 'names': 'gene',
                                   'logfoldchanges': 'avg_logFC',
                                   'pvals': 'p_val', 'pvals_adj': 'p_val_adj', 'scores': 'scores',
                                   'pct_nz_group': 'pct.1', 'pct_nz_reference': 'pct.2'}, inplace=True)
        markers[res] = list(markers_ta.to_dict(orient='index').values())
        for x in markers[res]:
            x['avg_logFC'] = np.round(x['avg_logFC'], 4)
    with open(outfile, 'w') as out_file:
        json.dump(markers, out_file, indent=4, sort_keys=False)
