rule t_cells_Guj8Kie5:
    input: h5ad=rules.integration.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Guj8Kie5/data.h5ad'
    params: res=config['res'], cells_to_filter=config['t_cells_Guj8Kie5_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Guj8Kie5/meta_data.csv'
    log: 'logs/scanpy/t_cells_Guj8Kie5.log'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_t_cells_Guj8Kie5:
    input: h5ad=rules.t_cells_Guj8Kie5.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Guj8Kie5/t_cells_Guj8Kie5")
    params: token='t_cells_Guj8Kie5', name='t_cells_Guj8Kie5',
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
           
rule t_cells_Xa5Aixae:
    input: h5ad=rules.t_cells_Guj8Kie5.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Xa5Aixae/data.h5ad'
    params: res=config['res'], cells_to_filter=config['t_cells_Xa5Aixae_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Xa5Aixae/meta_data.csv'
    log: 'logs/scanpy/t_cells_Xa5Aixae.log'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_t_cells_Xa5Aixae:
    input: h5ad=rules.t_cells_Xa5Aixae.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Xa5Aixae/t_cells_Xa5Aixae")
    params: token='t_cells_Xa5Aixae', name='t_cells_Xa5Aixae',
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
           
rule t_cells_Ux4ooGha:
    input: h5ad=rules.t_cells_Xa5Aixae.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Ux4ooGha/data.h5ad'
    params: res=config['res'], cells_to_filter=config['t_cells_Ux4ooGha_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Ux4ooGha/meta_data.csv'
    log: 'logs/scanpy/t_cells_Ux4ooGha.log'
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

rule conversion_t_cells_Ux4ooGha:
    input: h5ad=rules.t_cells_Ux4ooGha.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/t_cells/t_cells_Ux4ooGha/t_cells_Ux4ooGha")
    params: token='t_cells_Ux4ooGha', name='t_cells_Ux4ooGha',
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
           
rule cd4_t_cells_Ux4ooGha:
    input: h5ad=rules.t_cells_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_Ux4ooGha/data.h5ad'
    params: res=config['res'], cells_to_filter=config['cd4_t_cells_Ux4ooGha_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_Ux4ooGha/meta_data.csv'
    log: 'logs/scanpy/cd4_t_cells_Ux4ooGha.log'
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

rule conversion_cd4_t_cells_Ux4ooGha:
    input: h5ad=rules.cd4_t_cells_Ux4ooGha.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_Ux4ooGha/cd4_t_cells_Ux4ooGha")
    params: token='cd4_t_cells_Ux4ooGha', name='cd4_t_cells_Ux4ooGha',
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
           
rule cd8_dn_t_cells_Ux4ooGha:
    input: h5ad=rules.t_cells_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd8_dn_t_cells/cd8_dn_t_cells_Ux4ooGha/data.h5ad'
    params: res=config['res'], cells_to_filter=config['cd8_dn_t_cells_Ux4ooGha_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd8_dn_t_cells/cd8_dn_t_cells_Ux4ooGha/meta_data.csv'
    log: 'logs/scanpy/cd8_dn_t_cells_Ux4ooGha.log'
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

rule conversion_cd8_dn_t_cells_Ux4ooGha:
    input: h5ad=rules.cd8_dn_t_cells_Ux4ooGha.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd8_dn_t_cells/cd8_dn_t_cells_Ux4ooGha/cd8_dn_t_cells_Ux4ooGha")
    params: token='cd8_dn_t_cells_Ux4ooGha', name='cd8_dn_t_cells_Ux4ooGha',
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
           
rule cd4_t_cells_ungi0oaP:
    input: h5ad=rules.cd4_t_cells_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_ungi0oaP/data.h5ad'
    params: res=config['res'], cells_to_filter=config['cd4_t_cells_ungi0oaP_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_ungi0oaP/meta_data.csv', correct_meta='correct_meta_data.csv'
    log: 'logs/scanpy/cd4_t_cells_ungi0oaP.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --corr_meta {params.correct_meta} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_cd4_t_cells_ungi0oaP:
    input: h5ad=rules.cd4_t_cells_ungi0oaP.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_ungi0oaP/cd4_t_cells_ungi0oaP")
    params: token='cd4_t_cells_ungi0oaP', name='cd4_t_cells_ungi0oaP',
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
           
rule cd4_t_cells_uiSae8ee:
    input: h5ad=rules.cd4_t_cells_ungi0oaP.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_uiSae8ee/data.h5ad'
    params: res=config['res'], cells_to_filter=config['cd4_t_cells_uiSae8ee_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_uiSae8ee/meta_data.csv', correct_meta='correct_meta_data.csv'
    log: 'logs/scanpy/cd4_t_cells_uiSae8ee.log'
    resources: tmpdir='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/tmp'
    threads: 6
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --corr_meta {params.correct_meta} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_cd4_t_cells_uiSae8ee:
    input: h5ad=rules.cd4_t_cells_uiSae8ee.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd4_t_cells/cd4_t_cells_uiSae8ee/cd4_t_cells_uiSae8ee")
    params: token='cd4_t_cells_uiSae8ee', name='cd4_t_cells_uiSae8ee',
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
           
rule cd8_dn_t_cells_uiSae8ee:
    input: h5ad=rules.cd8_dn_t_cells_Ux4ooGha.output.h5ad
    output: h5ad='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd8_dn_t_cells/cd8_dn_t_cells_uiSae8ee/data.h5ad'
    params: res=config['res'], cells_to_filter=config['cd8_dn_t_cells_uiSae8ee_to_filter'], meta='/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd8_dn_t_cells/cd8_dn_t_cells_uiSae8ee/meta_data.csv', corr_meta='correct_meta_data.csv'
    log: 'logs/scanpy/cd8_dn_t_cells_uiSae8ee.log'
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

rule conversion_cd8_dn_t_cells_uiSae8ee:
    input: h5ad=rules.cd8_dn_t_cells_uiSae8ee.output.h5ad
    output: scn_dir=directory("/storage1/fs1/martyomov/Active/collaborations/greg_brian/greg_wu/csf_atlas/output/cd8_dn_t_cells/cd8_dn_t_cells_uiSae8ee/cd8_dn_t_cells_uiSae8ee")
    params: token='cd8_dn_t_cells_uiSae8ee', name='cd8_dn_t_cells_uiSae8ee',
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