rule reanalysis_Bib2vieg:
    input: h5ad=rules.integration.output.h5ad
    output: h5ad='output/csf_tabula_Bib2vieg/data.h5ad',
            meta='output/csf_tabula_Bib2vieg/meta_data.csv'
    params: res=config['res'], cells_to_filter=config['csf_tabula_Bib2vieg_to_filter']
    log: 'logs/scanpy/csf_tabula_Bib2vieg.log'
    threads: 12
    shell: """
           python workflow/scripts/reanalysis.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {output.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_filtered:
    input: h5ad=rules.reanalysis_Bib2vieg.output.h5ad
    output: scn_dir="output/csf_tabula_Bib2vieg/csf_tabula_Bib2vieg"
    params: token='csf_tabula_Bib2vieg', name='csf_tabula_Bib2vieg',
            specie=config['specie']
    shell: """
           python workflow/scripts/conversion.py \
           --adata {input.h5ad} \
           --token {params.token} \
           --name {params.name} \
           --out_dir {output.scn_dir} \
           --specie {params.specie}
           """
           
           
rule b_cells_uW0iG3oh:
    input: h5ad=rules.reanalysis_Bib2vieg.output.h5ad
    output: h5ad='output/b_cells_uW0iG3oh/data.h5ad'
    params: res=config['res'], cells_to_filter=config['b_cells_uW0iG3oh_to_filter'], meta='output/b_cells_uW0iG3oh/meta_data.csv'
    log: 'logs/scanpy/b_cells_uW0iG3oh.log'
    threads: 12
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {params.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_b_cells_uW0iG3oh:
    input: h5ad=rules.b_cells_uW0iG3oh.output.h5ad
    output: scn_dir=directory("output/b_cells_uW0iG3oh/b_cells_uW0iG3oh")
    params: token='b_cells_uW0iG3oh', name='b_cells_uW0iG3oh',
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
           
rule t_cells_uW0iG3oh:
    input: h5ad=rules.reanalysis_Bib2vieg.output.h5ad
    output: h5ad='output/t_cells/t_cells_uW0iG3oh/data.h5ad',
            meta='output/t_cells/t_cells_uW0iG3oh/meta_data.csv'
    params: res=config['res'], cells_to_filter=config['t_cells_uW0iG3oh_to_filter']
    log: 'logs/scanpy/t_cells_uW0iG3oh.log'
    threads: 12
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {output.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_t_cells_uW0iG3oh:
    input: h5ad=rules.t_cells_uW0iG3oh.output.h5ad
    output: scn_dir=directory("output/t_cells/t_cells_uW0iG3oh/t_cells_uW0iG3oh")
    params: token='t_cells_uW0iG3oh', name='t_cells_uW0iG3oh',
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
           
rule myeloid_cells_uW0iG3oh:
    input: h5ad=rules.reanalysis_Bib2vieg.output.h5ad
    output: h5ad='output/myeloid_cells/myeloid_cells_uW0iG3oh/data.h5ad',
            meta='output/myeloid_cells/myeloid_cells_uW0iG3oh/meta_data.csv'
    params: res=config['res'], cells_to_filter=config['myeloid_cells_uW0iG3oh_to_filter']
    log: 'logs/scanpy/myeloid_cells_uW0iG3oh.log'
    threads: 12
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {output.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """

rule conversion_myeloid_cells_uW0iG3oh:
    input: h5ad=rules.myeloid_cells_uW0iG3oh.output.h5ad
    output: scn_dir=directory("output/myeloid_cells/myeloid_cells_uW0iG3oh/myeloid_cells_uW0iG3oh")
    params: token='myeloid_cells_uW0iG3oh', name='myeloid_cells_uW0iG3oh',
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

           
rule nk_cells_uW0iG3oh:
    input: h5ad=rules.reanalysis_Bib2vieg.output.h5ad
    output: h5ad='output/nk_cells/nk_cells_uW0iG3oh/data.h5ad',
            meta='output/nk_cells/nk_cells_uW0iG3oh/meta_data.csv'
    params: res=config['res'], cells_to_filter=config['nk_cells_uW0iG3oh_to_filter']
    log: 'logs/scanpy/nk_cells_uW0iG3oh.log'
    threads: 12
    shell: """
           python workflow/scripts/reanalysis_by_study.py \
           --data {input.h5ad} \
           --out_file {output.h5ad} \
           --meta_out_file {output.meta} \
           --res {params.res} \
           --cells_to_filter {params.cells_to_filter} &> {log}
           """
           
rule conversion_nk_cells_uW0iG3oh:
    input: h5ad=rules.nk_cells_uW0iG3oh.output.h5ad
    output: scn_dir=directory("output/nk_cells/nk_cells_uW0iG3oh/nk_cells_uW0iG3oh")
    params: token='nk_cells_uW0iG3oh', name='nk_cells_uW0iG3oh',
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