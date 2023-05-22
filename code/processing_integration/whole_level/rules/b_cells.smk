rule b_cells_Guj8Kie5:
    input: h5ad=rules.integration.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Guj8Kie5/data.h5ad'
    params: res=config['res'], cells_to_filter=config['b_cells_Guj8Kie5_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Guj8Kie5/meta_data.csv'
    log: 'logs/scanpy/b_cells_Guj8Kie5.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_b_cells_Guj8Kie5:
    input: h5ad=rules.b_cells_Guj8Kie5.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Guj8Kie5/b_cells_Guj8Kie5")
    params: token='b_cells_Guj8Kie5', name='b_cells_Guj8Kie5',
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
           
           
rule b_cells_Xa5Aixae:
    input: h5ad=rules.b_cells_Guj8Kie5.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Xa5Aixae/data.h5ad'
    params: res=config['res'], cells_to_filter=config['b_cells_Xa5Aixae_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Xa5Aixae/meta_data.csv'
    log: 'logs/scanpy/b_cells_Xa5Aixae.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_b_cells_Xa5Aixae:
    input: h5ad=rules.b_cells_Xa5Aixae.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Xa5Aixae/b_cells_Xa5Aixae")
    params: token='b_cells_Xa5Aixae', name='b_cells_Xa5Aixae',
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
           
rule b_cells_Ux4ooGha:
    input: h5ad=rules.b_cells_Xa5Aixae.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Ux4ooGha/data.h5ad'
    params: res=config['res'], cells_to_filter=config['b_cells_Ux4ooGha_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Ux4ooGha/meta_data.csv'
    log: 'logs/scanpy/b_cells_Ux4ooGha.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_b_cells_Ux4ooGha:
    input: h5ad=rules.b_cells_Ux4ooGha.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_Ux4ooGha/b_cells_Ux4ooGha")
    params: token='b_cells_Ux4ooGha', name='b_cells_Ux4ooGha',
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

rule b_cells_j2gVns0:
    input: h5ad=rules.b_cells_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_j2gVns0/data.h5ad'
    params: res=config['res'], cells_to_filter=config['b_cells_j2gVns0_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_j2gVns0/meta_data.csv', corr_meta='correct_meta_data.csv'
    log: 'logs/scanpy/b_cells_j2gVns0.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --corr_meta {params.corr_meta} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_b_cells_j2gVns0:
    input: h5ad=rules.b_cells_j2gVns0.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/b_cells/b_cells_j2gVns0/b_cells_j2gVns0")
    params: token='b_cells_j2gVns0', name='b_cells_j2gVns0',
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