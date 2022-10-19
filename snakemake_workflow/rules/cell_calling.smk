samplesheet = samplesheets

include: "preprocessing.smk"
include: "vdjc.smk"


rule whitelist_cell_barcodes:
    input:
        "{base}/deduplicate_sample/{sample_uid}_consensus-pass_seq_ids_barcodes.tsv"
    output:
        "{base}/whitelisted/{sample_uid}_consensus-pass_seq_ids_barcodes_whitelisted.tsv"
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
        WHITELIST=config["whitelist"]
    log:
        "{base}/logs/{sample_uid}_whitelist_cbs.log"
    shell:
        """
        python {params.scripts}/whitelist_cbs_firstpass.py {input}\
        -whitelist {params.WHITELIST} \
        -outdir {wildcards.base}/whitelisted \
        2> {log}
        """

def aggregate_BC_maps(wildcards):
    samplelist = [
        "{}/whitelisted/{}_consensus-pass_seq_ids_barcodes_whitelisted.tsv".format(wildcards.base, bc)
        for bc in samplesheets[samplesheets.donor == wildcards.sample].index
    ]
    return samplelist

rule call_cells:
    input:
        BC_maps=aggregate_BC_maps, 
        VDJ_file="{base}/vseq_aligned_to_germline_polished/{sample}_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz",
    output:
        called_cells="{base}/cell_calls/{sample}_called_cells.tsv.gz",
        ambient_rna="{base}/cell_calls/{sample}_ambient_vdjs.tsv.gz",
    log:
        "{base}/logs/{sample}_call_cells.log"
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
        samplesheet=config["samplesheets"],
    shell:
        """
        python {params.scripts}/call_cells.py {input.VDJ_file} \
        -cb_umi_seqid_map {input.BC_maps} \
        -outdir {wildcards.base}/cell_calls \
        -samplesheet {params.samplesheet} \
        2> {log}
        """
rule combine_cells_with_vdj_annotations:
    input:
        vdj_info="{base}/vseq_aligned_to_germline_polished/{sample}_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz",
        cell_calls="{base}/cell_calls/{sample}_called_cells.tsv.gz"
    output:
        "{base}/cell_calls/{sample}_called_cells_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz"
    log:
        "{base}/logs/{sample}_combined_cell_call_and_annotation.log"
    params:
        scripts=config["scripts"]
    shell:
        """
        python {params.scripts}/combine_called_cell_with_vseq_data.py \
        -cell_calls {input.cell_calls} \
        -vdj_info {input.vdj_info} \
        -outdir {base}/cell_calls \
        2> {log}
        """
