rule get_cells:
    input: rda='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/data/{dataset}/object.RData'
    output: dir=directory('output/counts/{dataset}')
    params: name=lambda wildcards, input: wildcards.dataset,
            json=config['json'], annot=config['annot']
    threads: 1
    shell: """
           Rscript workflow/scripts/get_data.R \
           --rda {input.rda} \
           --json {params.json} \
           --out_dir {output.dir} \
           --name {params.name} \
           --annot {params.annot}
           """
           
           
rule integration:
    input: expand('output/counts/{dataset}', dataset=datasets.keys())
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/csf_tabula_agael6aC/data.h5ad',
            meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/csf_tabula_agael6aC/meta_data.csv'
    params: res=config['res'], atlas_json=config['atlas_json']
    log: 'logs/integration.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/csf_tabula_agael6aC/tmp'
    threads: 12
    shell: """
           python workflow/scripts/integration.py \
           --out_file {output.h5ad} \
           --meta_out_file {output.meta} \
           --res {params.res} \
           --config {params.atlas_json} 2> {log}
           """

rule conversion:
    input: h5ad=rules.integration.output.h5ad
    output: scn_dir=directory('/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/csf_tabula_agael6aC/csf_tabula_agael6aC/csf_tabula_agael6aC/csf_tabula_agael6aC/')
    params: token=config['token'], name=config['name'],
            specie=config['specie'], meta='correct_meta_data.csv'
    shell: """
           python workflow/scripts/conversion_with_correct_meta.py \
           --adata {input.h5ad} \
           --meta {params.meta} \
           --token {params.token} \
           --name {params.name} \
           --out_dir {output.scn_dir} \
           --specie {params.specie}
           """