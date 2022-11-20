import sys
import os

# -----------------------------
# Configuration

species = config["species"]

rule cellranger_count:
    output:
        "{base}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
        directory("{base}/per_sample/cellranger/{sample_uid}"),
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

rule run_decontx_per_sample:
    input:
        ancient("{base}/per_sample/cellranger/{sample_uid}"),
    output:
        "{base}/per_sample/decontX/{sample_uid}.h5ad",
    container:
        "docker://campbio/sctk_qc:1.7.6"
    params:
        scripts=config["scripts"],
        output_dir="{base}/decontX/",
        mode=config["decontX_mode"],
    resources:
        mem_mb="32000",
        partition="owners,quake",
        time="0-2",
    log:
        "{base}/logs/decontX/{sample_uid}_decontX.log",
    shell:
        "Rscript {params.scripts}/post_cellranger/decontX.R"
        " -i {input}"
        " -o {params.output_dir}"
        " -n {wildcards.sample_uid}"
        " -m {params.mode}"

rule add_sampleuid_column:
    input:
        "{base}/per_sample/decontX/{sample_uid}.h5ad",
    output:
        "{base}/per_sample/add_column/{sample_uid}.h5ad"
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    resources:
        mem_mb="16000",
        partition="owners,quake",
        time="0-2",
    log:
        "{base}/logs/add_column/{sample_uid}.log",
    script:        
        "../scripts/post_cellranger/add_sampleuid_column.py"

rule aggregate_h5ads:
    """ aggreate h5ads converted from cellranger """
    input:
        expand(
            "{base}/per_sample/add_column/{sample_uid}.h5ad",
            base=config["base"],
            sample_uid=sample_uids,
        ),
    output:
        "{base}/aggregated/aggr_gex.h5ad",
    log:
        "{base}/logs/aggregation.log",
    resources:
        mem_mb="210000",
        partition="quake",
        time="0-1",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=3,
        min_counts=10,
        filter_cells = False
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_with_scanpy.py"

rule preprocess_scanpy:
    """ performs some preprocessing and qc as well as celltypist labeling, adds samplesheet info to object"""
    input:
        "{base}/aggregated/aggr_gex.h5ad",
        config["samplesheets"],
    output:
        "{base}/analysis/scanpy/gex_object.h5ad",
        "{base}/analysis/scanpy/assayFile.csv",
        "{base}/analysis/scanpy/annotFile.csv",
        "{base}/analysis/scanpy/featureFile.csv",
    resources:
        mem_mb="64000",
        partition="owners,quake",
        time="0-2",
    log:
        "{base}/logs/scanpy/preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/scanpy_pp.py"
