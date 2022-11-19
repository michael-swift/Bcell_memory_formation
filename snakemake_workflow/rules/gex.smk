import sys
import os

# -----------------------------
# Configuration

species = config["species"]

rule cellranger_count:
    output:
        "{base}/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
        directory("{base}/cellranger/{sample_uid}"),
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

rule filter_h5_write_h5ad:
    """ converts txg outputs with minimal filtering, adds column called sample_uid"""
    input:
        "{base}/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
    output:
        "{base}/per_sample/h5ads/filtered_txg_{sample_uid}.h5ad",
    log:
        "{base}/logs/{sample_uid}_aggregation.log",
    resources:
        mem_mb="32000",
        partition="normal,owners,quake",
        time="0-1",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=5,
        min_counts=10,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/h5_to_h5ad.py"

rule aggregate_h5ads:
    """ aggreate h5ads converted from cellranger """
    input:
        expand(
            "{base}/per_sample/h5ads/filtered_txg_{sample_uid}.h5ad",
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
    log:
        "{base}/logs/scanpy/preprocess.log",
    resources:
        mem_mb="64000",
        partition="owners,quake",
        time="0-2",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/scanpy_pp.py"

rule run_decontx:
    input:
        "{base}/analysis/scanpy",
        "{base}/analysis/scanpy/assayFile.csv",
        "{base}/analysis/scanpy/annotationFile.csv",
        "{base}/analysis/scanpy/featureFile.csv",
    output:
        directory("{base}/decontX"),
        "{base}/decontX/gex_object.h5ad",
    container:
        "docker://campbio/sctk_qc:1.7.6"
    params:
        scripts=config["scripts"],
        output_dir="{base}/decontX/",
        mode=config["decontX_mode"],
        label="celltypist",
        name="gex_object"
    resources:
        mem_mb="128000",
        partition="owners,quake",
        time="0-2",
    log:
        "{base}/logs/decontX/decontX.log",
    shell:
        "Rscript {params.scripts}/post_cellranger/decontX_all.R"
        " -c {input[1]}"
        " -a {input[2]}"
        " -f {input[3]"
        " -o {output[0]}"
        " -l {params.label}"
        " -n {params.name}"
        " -m {params.mode}"

rule run_decontx_per_sample:
    input:
        ancient("{base}/cellranger/{sample_uid}"),
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





















"""
rule preprocess_scanvi:
    input:
        "{base}/scanpy/aggr_gex.h5ad",
        config["samplesheets"],
    output:
        "{base}/analysis/scanvi/gex_object.h5ad.gz",
        "{base}/analysis/scanvi/gex_object_bcells.h5ad.gz",
    log:
        "{base}/logs/scanvi/preprocess.log",
    resources:
        mem_mb="210000",
        disk_mb="120000",
        partition="quake",
        time="0-2",
    conda:
        config["workflow_dir"] + "/envs/scvi.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/scanvi_pp.py"

rule aggregate_decontX:
    "Aggregates decontX objects"
    input:
        expand(
            "{base}/decontX/{sample_uid}.h5ad",
            base=config["base"],
            sample_uid=sample_uids,
        ),
    output:
        "{base}/aggregated/decontX_aggr_gex.h5ad",
    log:
        "{base}/logs/aggregation.log",
    resources:
        mem_mb="210000",
        partition="quake",
        time="0-1",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    params:
        min_genes=300,
        min_counts=100,
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_with_scanpy.py"
"""
