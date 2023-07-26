smplesheet = samplesheets
include: "gex.smk"
include: "annotate.smk"
include: "cell_calling.smk"

rule integrate_gex_and_vdj_data:
    input:
        expand("{base_vdj}/all_vdj_cell_calls_IGH.tsv.gz", base_vdj = base['vdj']),
        expand("{base_gex}/annotate/all_adata.obs.tsv.gz", base_gex = base['gex'])
    output:
        expand("{base_vdj}/integrated_cell_calls.tsv.gz", base_vdj = base['vdj'])
    params:
        scripts=config["vdj_scripts"],
        base_vdj=config["base"]['vdj'],
        base_gex=config["base"]['gex'],
    log:
        expand("{base_vdj}/logs/integrate_gex_and_vdj_data.log", base_vdj = base['vdj']),
    resources:
        mem_meb="65000",
    conda:
        "../envs/scanpy.py",
    shell:
        "python {params.scripts}/integrate_gex_and_vdj.py "
        "-vdj {input.vdj} "
        "-gex {input.gex} "
        "-output {output} "
        "2> {log}"
