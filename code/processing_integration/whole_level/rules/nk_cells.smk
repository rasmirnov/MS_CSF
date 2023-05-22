rule nk_cells_Guj8Kie5:
    input: h5ad=rules.integration.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/nk_cells/nk_cells_Guj8Kie5/data.h5ad'
    params: res=config['res'], cells_to_filter=config['nk_cells_Guj8Kie5_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/nk_cells/nk_cells_Guj8Kie5/meta_data.csv'
    log: 'logs/scanpy/nk_cells_Guj8Kie5.log'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_nk_cells_Guj8Kie5:
    input: h5ad=rules.nk_cells_Guj8Kie5.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/nk_cells/nk_cells_Guj8Kie5/nk_cells_Guj8Kie5")
    params: token='nk_cells_Guj8Kie5', name='nk_cells_Guj8Kie5',
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
           
           
rule nk_cells_Xa5Aixae:
    input: h5ad=rules.nk_cells_Guj8Kie5.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/nk_cells/nk_cells_Xa5Aixae/data.h5ad'
    params: res=config['res'], cells_to_filter=config['nk_cells_Xa5Aixae_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/nk_cells/nk_cells_Xa5Aixae/meta_data.csv'
    log: 'logs/scanpy/nk_cells_Xa5Aixae.log'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_nk_cells_Xa5Aixae:
    input: h5ad=rules.nk_cells_Xa5Aixae.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/nk_cells/nk_cells_Xa5Aixae/nk_cells_Xa5Aixae")
    params: token='nk_cells_Xa5Aixae', name='nk_cells_Xa5Aixae',
            specie=config['specie'], meta='correct_meta_data.csv'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    shell: """
           python workflow/scripts/conversion_with_correct_meta.py \
           --adata {input.h5ad} \
           --meta {params.meta} \
           --token {params.token} \
           --name {params.name} \
           --out_dir {output.scn_dir} \
           --specie {params.specie}
           """