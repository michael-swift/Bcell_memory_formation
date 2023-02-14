def expected_cells_adj(wildcards, attempt):
    return 10000 / attempt


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
        "{base}/per_sample/cellranger/sample_stamps/{sample_uid}.done"
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        transcriptome=config["transcriptome"],
        local_cores = 20
    shell:
        "mkdir -p {base}/per_sample/cellranger/ && "
        "cd {base}/per_sample/cellranger/ && "
        "rm -f {wildcards.sample_uid}/_lock && "
        "{params.cell_ranger}/cellranger count"
        " --id={wildcards.sample_uid}"
        " --transcriptome={params.transcriptome}"
        " --fastqs={input}"
        " --localcores={params.local_cores}"
        " --sample={wildcards.sample_uid}"
        " --nosecondary"
        " --no-bam && touch {output}"

rule run_cellbender:
    input:rules.cellranger_count.output
    output:
        "{base}/per_sample/cellbender/{sample_uid}/background_removed.h5",
        "{base}/per_sample/cellbender/{sample_uid}/background_removed_cell_barcodes.csv",
    conda:
        "cellbender2"
    log:
        "{base}/logs/cellbender/{sample_uid}.log",
    resources:
        time="4:00:00",
        partition="gpu",
    params:
        expected_cells=lambda wildcards: samplesheet_lookup(wildcards.sample_uid, "expected_cells"),
        cellranger_count_file = lambda wildcards: "{base}/per_sample/{sample_uid}/outs/raw_feature_bc_matrix.h5".format(sample_uid = wildcards.sample_uid, base = config["base"])
    shell:
        "cellbender remove-background"
        " --input {params.cellranger_count_file}"
        " --output {output[0]}"
        " --expected-cells {params.expected_cells}"
        " --fpr 0.01"
        " --cuda"

rule combine_cb_cr:
    input:
        cr=rules.cellranger_count.output,
        cb="{base}/per_sample/cellbender/{sample_uid}/background_removed.h5"
    output:
        "{base}/per_sample/cellranger_cellbender/{sample_uid}/combined.h5ad"
    log:
        "{base}/logs/{sample_uid}/cb_cr_combined.log",
    resources:
        mem_mb="32000",
        partition="quake,owners",
        time="0-1",
    params:
        min_genes=200,
        min_counts=400,
        filter_cells=True,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/combined_cr_cb.py"

rule aggregate_h5ads:
    """ aggregate h5 from cellranger or h5ads scanpy """
    input:
        expand(
            "{base}/per_sample/cellranger_cellbender/{sample_uid}/combined.h5ad",
            base=config["base"],
            sample_uid=sample_uids,
        ),"{base}/downloads/CountAdded_PIP_global_object_for_cellxgene.h5ad"
    output:
        "{base}/aggregated/aggr_gex_raw.h5ad",
    log:
        "{base}/logs/aggregation.log",
    resources:
        mem_mb="210000",
        partition="quake,owners",
        time="0-1",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=300,
        min_counts=300,
        filter_cells=True,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_h5ads.py"

rule preprocess_scanpy:
    """ performs some preprocessing and qc as well as celltypist labeling, adds samplesheet info to object"""
    input:
        "{base}/aggregated/aggr_gex_raw.h5ad",
        config["samplesheets"],
    output:
        "{base}/analysis/scanpy/gex_object.h5ad.gz",
    resources:
        partition="quake,owners",
        mem_mb="210000",
    params:
        outdir="{base}/analysis/scanpy/h5ad_by_sample_uid",
    log:
        "{base}/logs/scanpy/preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/scanpy_pp.py"
