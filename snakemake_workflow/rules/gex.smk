import sys
import os

def expected_cells_adj(wildcards, attempt):
    return 10000 / attempt
# -----------------------------
# Configuration

species = config["species"]

rule cellranger_count:
    output:
        "{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        transcriptome=config["transcriptome"],
        fastq_dir=config["fastq_dir"],
    resources:
        mem_mb="210000",
        partition="quake",
        disk_mb="8000",
        time="2-0",
    threads: 20
    shell:
        "mkdir -p {params.base}/cellranger &&"
        " cd {params.base}/cellranger &&"
        " rm -rf {wildcards.sample_uid} && {params.cell_ranger}/cellranger count --id={wildcards.sample_uid} --transcriptome={params.transcriptome} --fastqs {params.fastq_dir} --sample={wildcards.sample_uid} --localcores=20 --nosecondary --no-bam > {log}"

rule run_cellbender:
    input:
        ancient("{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5")
    output:
        "{base}/per_sample/cellbender/{sample_uid}/background_removed.h5"
    conda:
        "cellbender2"
    log:
        "{base}/logs/cellbender/{sample_uid}.log",
    resources:
        time="4:00:00",
        partition="gpu",
        expected_cells=expected_cells_adj
    shell:
        "cellbender remove-background"
        " --input {input}"
        " --output {output}"
        " --expected-cells {resources.expected_cells}" # this should be in params but something about wildcards
        " --fpr 0.01"
        " --cuda"

rule aggregate_h5ads:
    """ aggreate h5 from cellranger or h5ads scanpy """
    input:
        expand(
                    "{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
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
        #"scanpy_latest"
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=10,
        min_counts=300,
        filter_cells = True
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_with_scanpy.py"

rule preprocess_scanpy:
    """ performs some preprocessing and qc as well as celltypist labeling, adds samplesheet info to object"""
    input:
        "{base}/aggregated/aggr_gex_raw.h5ad",
        config["samplesheets"],
    output:
        "{base}/analysis/scanpy/gex_object.h5ad",
    resources:
        partition="quake,owners",
        mem_mb="210000"
    params:
        outdir="{base}/analysis/scanpy/h5ad_by_sample_uids"
    log:
        "{base}/logs/scanpy/preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/scanpy_pp.py"

rule run_decontx:
    input:
        "{base}/analysis/scanpy/gex_object.h5ad",    
    output:
        "{base}/aggregated/decontX/gex_object.h5ad",
    container:
        "docker://campbio/sctk_qc:1.7.6"
    params:
        scripts=config["scripts"],
        output_dir="{base}/decontX/",
        mode=config["decontX_mode"],
    resources:
        mem_mb="210000",
        partition="quake",
        time="0-2",
    log:
        "{base}/logs/decontX/decontX.log",
    shell:
        "Rscript {params.scripts}/post_cellranger/decontX.R"
        " -i {input}"
        " -o {params.output_dir}"
        " -n "
        " -m {params.mode}"
