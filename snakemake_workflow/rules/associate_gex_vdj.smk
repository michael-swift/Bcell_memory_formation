smplesheet = samplesheets

include: "annotate.smk"
include: "cell_calling.smk"

rule integrate_gex_and_vdj_data:
    input:
        vdj="{params.base_vdj}/all_vdj_cell_calls_IGH.tsv.gz",
        gex="{params.base_gex}/annotate/adata.obs.tsv.gz",
    output:
        "{params.base_vdj}/integrated_cell_calls.tsv.gz",
    params:
        scripts=config["vdj_scripts"],
        base_vdj=config["base"]['vdj'],
        base_gex=config["base"]['gex'],
    log:
        "{params.base_vdj}/logs/integrate_gex_and_vdj_data.log",
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
