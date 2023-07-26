smplesheet = samplesheets
include: "gex.smk"
include: "annotate.smk"
include: "cell_calling.smk"

rule integrate_gex_and_vdj_data:
    input:
        vdj=expand("{base_vdj}/all_vdj_cell_calls_IGH.tsv.gz", base_vdj = base['vdj']),
        gex=expand("{base_gex}/annotate/all_adata.obs.tsv.gz", base_gex = base['gex'])
    output:
        expand("{base_vdj}/integrated_cell_calls.tsv.gz", base_vdj = base['vdj'])
    params:
        scripts=config["vdj_scripts"],
        info=config['sample_relationships'],
    log:
        expand("{base_vdj}/logs/integrate_gex_and_vdj_data.log", base_vdj = base['vdj']),
    resources:
        mem_mb="128000",
    conda:
        "../envs/pacbio.yaml",
    shell:
        "python {params.scripts}/integrate_gex_and_vdj.py "
        "-vdj {input.vdj} "
        "-gex {input.gex} "
        "-output {output} "
        "-info {params.info} "
        "2> {log}"
