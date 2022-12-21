samplesheet = samplesheets


def expected_cells_adj(wildcards, attempt):
    return 10000 / attempt


def wildcard_input(wildcards):
    return


# -----------------------------
# Configuration
#
species = config["species"]


rule cellranger_count:
    input:
        config["fastq_dirs"],
    output:
        "{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        transcriptome=config["transcriptome"],
    resources:
        partition="quake,owners",
        disk_mb="8000",
        time="1-12",
    threads: 20
    shell:
        "mkdir -p {base}/per_sample/cellranger && "
        "cd {base}/per_sample/cellranger && "
        "rm -rf {wildcards.sample_uid} && "
        "{params.cell_ranger}/cellranger count"
        " --id={wildcards.sample_uid}"
        " --transcriptome={params.transcriptome}"
        " --fastqs={input[0]}"
        " --sample={wildcards.sample_uid}"
        " --localcores=20"
        " --nosecondary"
        " --no-bam > {log}"

rule run_cellbender:
    input:"{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"
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
    shell:
        "cellbender remove-background"
        " --input {input}"
        " --output {output[0]}"
        " --expected-cells {params.expected_cells}"
        " --fpr 0.01"
        " --cuda > {log}"

rule convert_h5_to_h5ad:
    input:
        cr="{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
        cb="{base}/per_sample/cellbender/{sample_uid}/background_removed.h5"
    output:
        "{base}/per_sample/cellranger/{sample_uid}/TenX.h5ad",
        "{base}/per_sample/cellbender/{sample_uid}/CellBender.h5ad"
    log:
        "{base}/logs/{sample_uid}/cb_cr_h5ads.log",
    resources:
        mem_mb="32000",
        partition="quake,owners",
        time="0-1",
    conda:
        #"scanpy_latest"
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=100,
        min_counts=300,
        filter_cells=True,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/h5_to_h5ad.py > {log}"

rule combine_cellbender_cellranger:
    input:
        "{base}/per_sample/cellranger/{sample_uid}/TenX.h5ad",
        "{base}/per_sample/cellbender/{sample_uid}/CellBender.h5ad"
    output:
        "{base}/per_sample/cb_and_cr/{sample_uid}/combined.h5ad"
    log:
        "{base}/logs/{sample_uid}/cb_cr_h5ads.log",
    resources:
        mem_mb="32000",
        partition="quake,owners",
        time="0-1",
    conda:
        #"scanpy_latest"
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=100,
        min_counts=300,
        filter_cells=True,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/h5_to_h5ad.py > {log}"


rule aggregate_h5ads:
    """ aggregate h5 from cellranger or h5ads scanpy """
    input:
        expand(
            "{base}/per_sample/cb_and_cr/{sample_uid}/combined.h5ad",
            base=config["base"],
            sample_uid=sample_uids,
        ),
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
        min_genes=10,
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
        outdir="{base}/analysis/scanpy/h5ad_by_sample_uids",
    log:
        "{base}/logs/scanpy/preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/scanpy_pp.py"
