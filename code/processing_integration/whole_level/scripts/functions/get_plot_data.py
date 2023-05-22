import re
import json
import pandas as pd
from anndata import AnnData


def get_plot_data(adata: AnnData, outfile: str, max_reduction_dims=5, n_levels=50):
    plot_data = {'fields': {}, 'data': {}, 'annotations': {}}
    for idx, r in enumerate(adata.obsm):
        if idx:
            new_df = pd.DataFrame(adata.obsm[r]).iloc[:, :min(len(adata.obsm[r][0]), max_reduction_dims + 1)]
            new_df.columns = [f"{re.sub('X_', '', r)}_{NUM + 1}".upper() for NUM in range(0, len(new_df.columns))]
            df = pd.concat([df.reset_index(drop=True), new_df], axis=1)
        else:
            df = pd.DataFrame(adata.obsm[r]).iloc[:, :min(len(adata.obsm[r][0]), max_reduction_dims + 1)]
            df.columns = [f"{re.sub('X_', '', r)}_{NUM + 1}".upper() for NUM in range(0, len(df.columns))]
    df.index = adata.obs.index
    info = pd.concat([adata.obs, df], axis=1)
    info['_row'] = info.index
    info_dict = info.to_dict(orient='index')
    plot_data['data'] = list(info_dict.values())
    categorical_data = info.select_dtypes(include='category')
    for k in info.columns:
        if k in df.columns:
            plot_data['fields'][k] = {"type": "numeric", "range": [min(info[k]), max(info[k])]}
        else:
            if k in categorical_data.columns:
                if len(list(adata.obs[k].dtype.categories)) < n_levels:
                    plot_data['fields'][k] = {"type": 'factor',
                                              "levels": list(adata.obs[k].dtype.categories)}
                    centers: dict = pd.DataFrame.from_dict(info_dict, orient = 'index').groupby(k).median()[
                        ['UMAP_1', 'UMAP_2']].to_dict(orient='index')
                    plot_data['annotations'][f'umap_{k}_centers'] = {'type': 'text',
                                                                     'value': k,
                                                                     'coords': ['UMAP_1', 'UMAP_2'],
                                                                     'data': list()}
                    for col, v in centers.items():
                        center = v
                        center[k] = str(col)
                        center['Text'] = col
                        plot_data['annotations'][f'umap_{k}_centers']['data'].append(center)

                else:
                    raise ValueError(f'Categorical variable {k} has more than {n_levels} values, make sure you want to convert it')
            else:
                if k != '_row' and k not in plot_data['fields'].keys():
                    plot_data['fields'][k] = {"type": "numeric", "range": [min(info[k]), max(info[k])]}

    with open(outfile, 'w') as out_file:
        json.dump(plot_data, out_file, indent=4, sort_keys=False)
