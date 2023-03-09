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
        "{base_gex}/per_sample/cellranger/sample_stamps/{sample_uid}.done","{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"

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

rule run_cellbender:
    input:ancient("{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5")
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
    shell:
        "cellbender remove-background"
        " --input {input}"
        " --output {output[0]}"
        " --expected-cells {params.expected_cells}"
        " --fpr 0.01"
        " --cuda"

rule combine_cb_cr:
    input:
        cr=ancient("{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"),
        cb=ancient("{base_gex}/per_sample/cellbender/{sample_uid}/background_removed.h5"),
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
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/combined_cr_cb.py"

rule aggregate_h5ads:
    """ aggregate h5 from cellranger or h5ads scanpy """
    input:
        expand(
            "{base_gex}/per_sample/cellranger_cellbender/{sample_uid}/combined.h5ad.gz",
            base_gex=config["base"]["gex"],
            sample_uid=sample_uids,
        ),
        "{base_gex}/downloads/TissueAdded_CountAdded_PIP_global_object_for_cellxgene.h5ad.gz",
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

rule annotate_by_tissue:
    """ performs preprocessing and qc as well as celltypist labeling, cell cycle assignment, doublet calling, adds samplesheet info to object"""
    input:
        "{base_gex}/aggregated/aggr_gex_raw.h5ad.gz",
        config["samplesheets"],
    output:
        tissue_dir=directory("{base_gex}/annotate/tissue_objs/{tissue}/"),
        tissue_obj="{base_gex}/annotate/tissue_objs/{tissue}/annotated_processed.h5ad.gz",
        labels="{base_gex}/annotate/tissue_objs/{tissue}/gex_labels.tab",
    resources:
        partition="quake,owners",
    params:
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"]
        ),
    log:
        "{base_gex}/logs/annotate/{tissue}/annotate_by_tissue.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/annotate.py"

rule aggregate_scanpy:
    input:
        h5ads=expand(
            "{base_gex}/annotate/tissue_objs/{tissue}/annotated_processed.h5ad.gz",
            base_gex=base['gex'],
            tissue=tissues,
        ),
        obs_dfs=expand(
            "{base_gex}/annotate/tissue_objs/{tissue}/gex_labels.tab",
            base_gex=base['gex'],
            tissue=tissues,
        ),
    output:
        full_h5ad="{base_gex}/annotate/gex_object.h5ad.gz",
        gex_labels="{base_gex}/annotate/gex_labels.tab.gz",
    resources:
        partition="quake,owners",
    params:
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"]
        ),
    log:
        "{base_gex}/logs/annotate/obj_preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_annotated.py"
