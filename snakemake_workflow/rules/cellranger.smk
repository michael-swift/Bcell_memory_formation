import sys
import os

# -----------------------------
# Configuration

species = config["species"]

rule cellranger_count:
    input:
    output:
        "{base}/10X/{sample_uid}/outs/raw_feature_bc_matrix.h5",
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        transcriptome=config["transcriptome"],
        fastq_dir=config['fastq_dir']
    resources: mem_mb="210000",partition="quake",disk_mb="8000",time="2-0"
    threads:20
    shell:
        "mkdir -p {params.base}/10X &&"
        " cd {params.base}/10X &&"
        " rm -rf {wildcards.sample_uid} && {params.cell_ranger}/cellranger count --id={wildcards.sample_uid} --transcriptome={params.transcriptome} --fastqs {params.fastq_dir} --sample={wildcards.sample_uid} --localcores=20 --nosecondary --no-bam"

# use scanpy to convert molecule.h5 file
rule h5_to_h5ad:
    """ Uses scanpy to convert h5 to h5ad (anndata object) """
    input:"{base}/10X/{sample_uid}/outs/raw_feature_bc_matrix.h5"
    output:
        "{base}/post_10X/{sample_uid}/raw_gex.h5ad"
    log:
        "{base}/logs/{sample_uid}/h5_to_h5ad.log"
    resources: mem_mb="128000",partition="quake",time="0-1"
    conda:config["workflow_dir"]+"/envs/scanpy.yaml"
    script:config["workflow_dir"]+"/scripts/h5_to_h5ad.py"

rule aggregate_gex:
    """ Aggregates all the genes expression data from 10X with very lenient filtering """ 
    input:expand("{base}/10X/{sample_uid}/outs/raw_feature_bc_matrix.h5", base=config["base"], sample_uid=sample_uids),
    output:
        "{base}/post_10X/aggr_gex_raw.h5ad"
    log:
        "{base}/logs/aggregation.log"
    resources: mem_mb="128000",partition="quake",time="0-1"
    conda:config["workflow_dir"]+"/envs/scanpy.yaml"
    params:min_genes="30", min_counts="100"
    script:config["workflow_dir"]+"/scripts/aggregate_with_scanpy.py"

rule preprocess_gex:
    input:"{base}/post_10X/aggr_gex_raw.h5ad", config["samplesheets"]
    output:
        "{base}/post_10X/processed/gex_object.h5ad.gz", "{base}/post_10X/processed/gex_object_bcells.h5ad.gz", "{base}/post_10X/processed/gex_object_bcell_subsampled.h5ad.gz"
    log:
        "{base}/logs/aggregation.log"
    resources: mem_mb="128000",partition="quake",time="0-1"
    conda:config["workflow_dir"]+"/envs/scanpy.yaml"
    script:config["workflow_dir"]+"/scripts/preprocess_GEX.py"

rule run_decontx:
    input:
        "{base}/post_10X/{sample_uid}/raw_gex.h5ad"
    output:
        "{base}/post_10X/{sample_uid}/decontX.h5ad"
    container:
        "docker://campbio/sctk_qc:1.7.6"
    params:
        scripts=config["scripts"],
    log:
        "{base}/logs/{sample_uid}_decontX.log"
    shell:
        "Rscript decontX.R {params.scripts}/decontX.R"
        " -i input[0]"
        " -o output[0]"
        " -n {sample_uid}"
        "2> {log}"
