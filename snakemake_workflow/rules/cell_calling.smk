samplesheet = samplesheets

include: "preprocessing.smk"
include: "vdjc.smk"


rule call_cells:
    input:
        VDJ_file=rules.realigned_to_polished_germline.output,
    output:
        called_cells="{base}/aggregated/cell_calls/{donor}_called_cells.tsv.gz",
        ambient_rna="{base}/aggregated/cell_calls/{donor}_ambient_vdjs.tsv.gz",
    log:
        "{base}/logs/{donor}_call_cells.log"
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
        whitelist=config["whitelist"],
    shell:
        
        "python {params.scripts}/call_cells.py {input.VDJ_file} "
        "-whitelist {params.whitelist} "
        "-outdir {wildcards.base}/aggregated/cell_calls "
        "-outname {wildcards.donor} "
        "2> {log}"
        
rule combine_cells_with_vdj_annotations:
    input:
        vdj_info=rules.realigned_to_polished_germline.output,
        cell_calls=rules.call_cells.output.called_cells,
    output:
        "{base}/aggregated/cell_calls/{donor}_called_cells_vdj_annotated.tsv.gz"
    log:
        "{base}/logs/{donor}_combined_cell_call_and_annotation.log"
    params:
        scripts=config["scripts"]
    shell:
        "python {params.scripts}/combine_called_cell_with_vseq_data.py "
        "-cell_calls {input.cell_calls} "
        "-vdj_info {input.vdj_info} "
        "-outdir {wildcards.basebase}/aggregated/cell_calls "
        "-outname {wildcards.donor}_called_cells_vdj_annotated"
        "2> {log}"