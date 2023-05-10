def wildcard_input(wildcards):
    return


# -----------------------------
# Configuration
#
species = config["species"]


rule cellranger_count:
    # snakemake and cellranger don't play well together because they both want control of the directory
    # touching .done file  to work around
    # unfortunately I want to refer to CellRanger Generated Files later in the workflow
    input:
        config["gex_fastq_dirs"],
    output:
        "{base_gex}/per_sample/cellranger/sample_stamps/{sample_uid}.done",
    # "{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
    params:
        name="count_cellranger",
        base=config["base"]["gex"],
        cell_ranger=config["cell_ranger"],
        transcriptome=config["transcriptome"],
    shell:
        "mkdir -p {wildcards.base_gex}/per_sample/cellranger/ && "
        "cd {wildcards.base_gex}/per_sample/cellranger/ && "
        "rm -f {wildcards.sample_uid}/_lock && "
        "{params.cell_ranger}/cellranger count"
        " --id={wildcards.sample_uid}"
        " --transcriptome={params.transcriptome}"
        " --fastqs={input}"
        " --sample={wildcards.sample_uid}"
        " --nosecondary"
        " --no-bam && touch {output}"


rule copy_cellranger:
    input:
        rules.cellranger_count.output,
    output:
        temp("{base_gex}/per_sample/cellranger_cp/{sample_uid}/outs/raw_feature_bc_matrix.h5")
    shell:
        "cp {wildcards.base_gex}/per_sample/cellranger/{wildcards.sample_uid}/outs/raw_feature_bc_matrix.h5 {output}"

rule run_cellbender:
    input:
            "{base_gex}/per_sample/cellranger_cp/{sample_uid}/outs/raw_feature_bc_matrix.h5",
    output:
        "{base_gex}/per_sample/cellbender/{sample_uid}/background_removed.h5",
        "{base_gex}/per_sample/cellbender/{sample_uid}/background_removed_cell_barcodes.csv",
    conda:
        "cellbender2"
    log:
        "{base_gex}/logs/cellbender/{sample_uid}.log",
    resources:
        time="4:00:00",
        partition="gpu",
    params:
        expected_cells=lambda wildcards: samplesheet_lookup(
            wildcards.sample_uid, "expected_cells"
        ),
        total_droplets=lambda wildcards: samplesheet_lookup(
            wildcards.sample_uid, "expected_cells"
        )
        * 3,
    shell:
        "cellbender remove-background"
        " --input {input}"
        " --output {output[0]}"
        " --total-droplets-included {params.total_droplets}"
        " --expected-cells {params.expected_cells}"
        " --fpr 0.01"
        " --cuda"


rule combine_cb_cr:
    input:
        cr=ancient(
            "{base_gex}/per_sample/cellranger_cp/{sample_uid}/outs/raw_feature_bc_matrix.h5"
        ),
        cb=ancient(
            "{base_gex}/per_sample/cellbender/{sample_uid}/background_removed.h5"
        ),
    output:
        "{base_gex}/per_sample/cellranger_cellbender/{sample_uid}/combined.h5ad.gz",
    log:
        "{base_gex}/logs/{sample_uid}/cb_cr_combined.log",
    resources:
        mem_mb="32000",
        partition="quake,owners",
        time="0-1",
    params:
        min_genes=1,
        min_counts=1,
        filter_cells=True,
        samplesheets=config["samplesheets"],
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/combine_cr_cb.py"


rule aggregate_h5ads:
    """ aggregate h5 from cellranger or h5ads scanpy """
    input:
        expand(
            "{base_gex}/per_sample/cellranger_cellbender/{sample_uid}/combined.h5ad.gz",
            base_gex=config["base"]["gex"],
            sample_uid=sample_uids,
        ),
        "{base_gex}/downloads/mod_PIP_global_object_for_cellxgene.h5ad.gz",
    output:
        "{base_gex}/aggregated/aggr_gex_raw.h5ad.gz",
    log:
        "{base_gex}/logs/aggregation.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=10,
        min_counts=50,
        filter_cells=True,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_h5ads.py"
