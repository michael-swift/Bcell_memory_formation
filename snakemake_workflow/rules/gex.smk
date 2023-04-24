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
        total_droplets=lambda wildcards: samplesheet_lookup(
            wildcards.sample_uid, "expected_cells"
        ) * 3
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
            "{base_gex}/per_sample/cellranger/{sample_uid}/outs/raw_feature_bc_matrix.h5"
        ),
        cb=ancient("{base_gex}/per_sample/cellbender/{sample_uid}/background_removed.h5")
        ,
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


rule pre_annotate_by_tissue:
    """ performs preprocessing and qc as well as celltypist labeling, cell cycle assignment, doublet calling, adds samplesheet info to object"""
    input:
        "{base_gex}/aggregated/aggr_gex_raw.h5ad.gz",
    output:
        fig_dir=directory("{base_gex}/pre_annotate/figures_{tissue}/"),
        tissue_obj="{base_gex}/pre_annotate/tissue_objs/{tissue}/pre_annotated_processed.h5ad.gz",
    resources:
        partition="quake,owners",
    params:
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"]
        ),
    log:
        "{base_gex}/logs/pre_annotate/{tissue}/pre_annotate_by_tissue.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/pre_annotate.py"


rule aggregate_pre_annotated:
    input:
        h5ads=expand(
            "{base_gex}/pre_annotate/tissue_objs/{tissue}/pre_annotated_processed.h5ad.gz",
            base_gex=base["gex"],
            tissue=tissues,
        ),
    output:
        full_h5ad="{base_gex}/pre_annotate/gex_object.h5ad.gz",
        adata_obs="{base_gex}/pre_annotate/adata.obs.tab.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/pre_annotate/obj_preprocess.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/post_cellranger/aggregate_pre_annotated.py"

rule remove_nonb:
    input:
        full_h5ad="{base_gex}/pre_annotate/gex_object.h5ad.gz",
    output:
        bcells="{base_gex}/pre_annotate/bcells.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/pre_annotate/removenonb_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
    run:
        import scanpy as sc
        adata = sc.read_h5ad(str(input))
        bcells = adata[adata.obs.probable_hq_single_b_cell.astype(str) == "True"]
        bcells.obs.loc[:,"atlas"] = adata.obs.donor.str.contains("TBd")
        bool_mapper = {True: "Tabula_Bursa", False: "TICA"}
        bcells.obs.loc[:, "atlas"] = bcells.obs["atlas"].map(bool_mapper)
        bcells.write_h5ad(output.bcells, compression="gzip")

rule scvi_bcells:
    input:
        bcell_h5ad="{base_gex}/pre_annotate/bcells.h5ad.gz",
    output:
        model=directory("{base_gex}/pre_annotate/scvi/model/covariates/"),
        adata="{base_gex}/pre_annotate/scvi/bcells.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/pre_annotate/scvi.obj_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"],
        ), subsample="False"
    conda:
        config["workflow_dir"] + "/envs/scvi.yaml"
    shell:
        "python {params.scripts}build_scvi.py {input} -output_file {output.adata} -output_model {output.model} -covariate_genes {params.cell_cycle_genes} -subsample {params.subsample}"

rule subset_bcells:
    input:
        adata="{base_gex}/pre_annotate/scvi/bcells.h5ad.gz"
    output:
        asc="{base_gex}/annotate/asc.h5ad.gz",
        mb="{base_gex}/annotate/mb.h5ad.gz",
        nb_other="{base_gex}/annotate/nb_other.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/pre_annotate/subset_bcells_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
    run:
        import scanpy as sc
        adata = sc.read_h5ad(str(input))
        adata = adata[adata.obs.Immune_All_Low_conf_score > 0.95]
        mb = adata[adata.obs.Immune_All_Low_predicted_labels.str.contains("Memory|Age")]
        mb.write_h5ad(output.mb, compression="gzip")
        asc = adata[adata.obs.Immune_All_Low_predicted_labels.str.contains("Plasma")]
        asc.write_h5ad(output.asc, compression="gzip")
        nb_other = adata[~adata.obs.Immune_All_Low_predicted_labels.str.contains("Plasma|Memory|Age")]
        nb_other.write_h5ad(output.nb_other, compression="gzip")
        
