rule myeloid_cells_Guj8Kie5:
    input: h5ad=rules.integration.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Guj8Kie5/data.h5ad'
    params: res=config['res'], cells_to_filter=config['myeloid_cells_Guj8Kie5_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Guj8Kie5/meta_data.csv'
    log: 'logs/scanpy/myeloid_cells_Guj8Kie5.log'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_myeloid_cells_Guj8Kie5:
    input: h5ad=rules.myeloid_cells_Guj8Kie5.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Guj8Kie5/myeloid_cells_Guj8Kie5")
    params: token='myeloid_cells_Guj8Kie5', name='myeloid_cells_Guj8Kie5',
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
           
           
rule myeloid_cells_Xa5Aixae:
    input: h5ad=rules.myeloid_cells_Guj8Kie5.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Xa5Aixae/data.h5ad'
    params: res=config['res'], cells_to_filter=config['myeloid_cells_Xa5Aixae_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Xa5Aixae/meta_data.csv'
    log: 'logs/scanpy/myeloid_cells_Xa5Aixae.log'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_myeloid_cells_Xa5Aixae:
    input: h5ad=rules.myeloid_cells_Xa5Aixae.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Xa5Aixae/myeloid_cells_Xa5Aixae")
    params: token='myeloid_cells_Xa5Aixae', name='myeloid_cells_Xa5Aixae',
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
           
rule myeloid_cells_Ux4ooGha:
    input: h5ad=rules.myeloid_cells_Xa5Aixae.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Ux4ooGha/data.h5ad'
    params: res=config['res'], cells_to_filter=config['myeloid_cells_Ux4ooGha_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Ux4ooGha/meta_data.csv'
    log: 'logs/scanpy/myeloid_cells_Ux4ooGha.log'
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

rule conversion_myeloid_cells_Ux4ooGha:
    input: h5ad=rules.myeloid_cells_Ux4ooGha.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/myeloid_cells/myeloid_cells_Ux4ooGha/myeloid_cells_Ux4ooGha")
    params: token='myeloid_cells_Ux4ooGha', name='myeloid_cells_Ux4ooGha',
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


rule microglia_Ux4ooGha:
    input: h5ad=rules.myeloid_cells_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_Ux4ooGha/data.h5ad'
    params: res=config['res'], cells_to_filter=config['microglia_Ux4ooGha_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_Ux4ooGha/meta_data.csv'
    log: 'logs/scanpy/microglia_Ux4ooGha.log'
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

rule conversion_microglia_Ux4ooGha:
    input: h5ad=rules.microglia_Ux4ooGha.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_Ux4ooGha/microglia_Ux4ooGha")
    params: token='microglia_Ux4ooGha', name='microglia_Ux4ooGha',
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
rule microglia_ungi0oaP:
    input: h5ad=rules.microglia_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_ungi0oaP/data.h5ad'
    params: res=config['res'], cells_to_filter=config['microglia_ungi0oaP_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_ungi0oaP/meta_data.csv', corr_meta='correct_meta_data.csv'
    log: 'logs/scanpy/microglia_ungi0oaP.log'
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

rule conversion_microglia_ungi0oaP:
    input: h5ad=rules.microglia_ungi0oaP.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_ungi0oaP/microglia_ungi0oaP")
    params: token='microglia_ungi0oaP', name='microglia_ungi0oaP',
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

rule microglia_ein8rohR:
    input: h5ad=rules.microglia_ungi0oaP.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_ein8rohR/data.h5ad'
    params: res=config['res'], cells_to_filter=config['microglia_ein8rohR_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_ein8rohR/meta_data.csv', corr_meta='correct_meta_data.csv'
    log: 'logs/scanpy/microglia_ein8rohR.log'
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

rule conversion_microglia_ein8rohR:
    input: h5ad=rules.microglia_ein8rohR.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_ein8rohR/microglia_ein8rohR")
    params: token='microglia_ein8rohR', name='microglia_ein8rohR',
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
           
rule microglia_joh5Iep2:
    input: h5ad=rules.microglia_ungi0oaP.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_joh5Iep2/data.h5ad'
    params: res=config['res'], cells_to_filter=config['microglia_joh5Iep2_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_joh5Iep2/meta_data.csv', corr_meta='correct_meta_data.csv'
    log: 'logs/scanpy/microglia_joh5Iep2.log'
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

rule conversion_microglia_joh5Iep2:
    input: h5ad=rules.microglia_joh5Iep2.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/microglia/microglia_joh5Iep2/microglia_joh5Iep2")
    params: token='microglia_joh5Iep2', name='microglia_joh5Iep2',
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
