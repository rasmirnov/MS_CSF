rule integration:
    output: rda="out/csf_tabula_ris_samples_ahk0buF0/object.RData"
    log: "logs/seurat/csf_tabula_ris_samples_ahk0buF0.log"
    params: out_dir='out/csf_tabula_ris_samples_ahk0buF0',
            sample_id='csf_tabula_ris_samples_ahk0buF0',
            js_f='annotations/integration.json'
    threads: 16
    shell: """
           Rscript workflow/scripts/integration.R \
           --js_f {params.js_f} \
           --out_dir {params.out_dir} \
           --sample_id {params.sample_id} 2> {log}
           """
           
rule integration_2:
    input: rda=rules.integration.output.rda
    output: rda="out/csf_tabula_ris_samples_Ierijuu6/object.RData"
    log: "logs/seurat/csf_tabula_ris_samples_Ierijuu6.log"
    params: out_dir='out/csf_tabula_ris_samples_Ierijuu6',
            sample_id='csf_tabula_ris_samples_Ierijuu6',
            js_f='annotations/integration_2.json'
    threads: 16
    shell: """
           Rscript workflow/scripts/reanalysis.R \
           --js_f {params.js_f} \
           --out_dir {params.out_dir} \
           --sample_id {params.sample_id} 2> {log}
           """