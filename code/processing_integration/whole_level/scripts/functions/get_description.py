import json

from anndata import AnnData


def get_description(adata: AnnData, token: str,
                    name: str, specie: str,
                    outfile: str, link: str, description: str,
                    public: bool, curated: bool):
    description = {'token': token,
                   'name': name,
                   'description': description,
                   'link': link,
                   'species': specie,
                   'cells': adata.obs.shape[0],
                   'public': public,
                   'curated': curated}
    with open(outfile, 'w') as out_file:
        json.dump(description, out_file, indent=4, sort_keys=False)
