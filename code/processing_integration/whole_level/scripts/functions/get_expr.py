import json

from anndata import AnnData

def get_expr(adata: AnnData, outfile: str):
    exp_data = {'features': list(adata.var.index),
     'barcodes': list(adata.obs.index),
     'totalCounts': list(adata.obs.nCount_RNA.astype(int)),
     'expType': 'counts'
    }
    with open(outfile, 'w') as out_file:
        json.dump(exp_data, out_file, indent=4, sort_keys=False)
