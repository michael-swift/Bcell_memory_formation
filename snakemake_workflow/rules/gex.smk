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
        "{base_gex}/per_sample/cellranger/sample_stamps/{sample_uid}.done",
        "{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5",
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
    input:
        ancient(
            "{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"
        ),
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
        cr=ancient(
            "{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"
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


rule annotate_by_tissue:
    """ performs preprocessing and qc as well as celltypist labeling, cell cycle assignment, doublet calling, adds samplesheet info to object"""
    input:
        "{base_gex}/aggregated/aggr_gex_raw.h5ad.gz",
    output:
        fig_dir=directory("{base_gex}/annotate/figures_{tissue}/"),
        tissue_obj="{base_gex}/annotate/tissue_objs/{tissue}/annotated_processed.h5ad.gz",
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


rule aggregate_annotated:
    input:
        h5ads=expand(
            "{base_gex}/annotate/tissue_objs/{tissue}/annotated_processed.h5ad.gz",
            base_gex=base["gex"],
            tissue=tissues,
        ),
    output:
        full_h5ad="{base_gex}/annotate/gex_object.h5ad.gz",
        adata_obs="{base_gex}/annotate/adata.obs.tab.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/obj_preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_annotated.py"


rule scvi_global:
    input:
        full_h5ad="{base_gex}/annotate/gex_object.h5ad.gz",
    output:
        model=directory("{base_gex}/annotate/scvi/model/covariates/"),
        adata="{base_gex}/annotate/scvi/gex_object.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/scvi.obj_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"],
        ), subsample="False"
    conda:
        config["workflow_dir"] + "/envs/scvi.yaml"
    shell:
        "python {params.scripts}build_scvi.py {input} -output_file {output.adata} -output_model {output.model} -covariate_genes {params.cell_cycle_genes} -subsample {params.subsample}"

rule bcell_subset:
    input:
        "{base_gex}/annotate/scvi/gex_object.h5ad.gz",
    output:
        bcells="{base_gex}/annotate/scvi/bcells.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/scvi.obj_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
    run:
        import scanpy as sc
        adata = sc.read_h5ad(str(input))
        bcells = adata[adata.obs.probable_hq_single_b_cell.astype(str) == "True"]
        sc.pp.neighbors(bcells, use_rep="X_scVI_cont", n_neighbors=20)
        sc.tl.umap(bcells, min_dist=0.3)
        sc.tl.leiden(bcells, key_added="leiden_scVI_cont", resolution=0.8)
        bcells.write_h5ad(output.bcells, compression="gzip")
