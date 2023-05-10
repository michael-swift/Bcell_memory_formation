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


rule scvi_all_cells:
    input:
        full_h5ad="{base_gex}/annotate/gex_object.h5ad.gz",
    output:
        adata="{base_gex}/annotate/scvi/all_cells.h5ad.gz",
    log:
        "{base_gex}/logs/annotate/scvi.obj_all.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"],
        ),
        subsample="False",
    conda:
        "scvi-2023"
    shell:
        "python {params.scripts}train_scvi_all.py {input} -output_file {output.adata} -covariate_genes {params.cell_cycle_genes} -subsample {params.subsample}"


rule remove_nonb:
    input:
        full_h5ad="{base_gex}/annotate/scvi/all_cells.h5ad.gz",
    output:
        bcells=temp("{base_gex}/annotate/bcells.h5ad.gz"),
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/removenonb_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
    run:
        import scanpy as sc

        adata = sc.read_h5ad(str(input))
        bcells = adata[adata.obs.probable_hq_single_b_cell.astype(str) == "True"]
        bcells.obs.loc[:, "atlas"] = adata.obs.donor.str.contains("TBd")
        bool_mapper = {True: "Tabula_Bursa", False: "TICA"}
        bcells.obs.loc[:, "atlas"] = bcells.obs["atlas"].map(bool_mapper)
        bcells.write_h5ad(output.bcells, compression="gzip")


rule pseudo_make_vdj:
    "need to connect to vdjc pipeline"
    output:
        "{base_gex}/vdjc/integrated_no_contaminants_with_all_ambient_flags.tsv.gz",
    shell:
        "cat {output}"


rule scvi_bcells:
    input:
        bcell_h5ad="{base_gex}/annotate/bcells.h5ad.gz",
    output:
        model=directory("{base_gex}/annotate/scvi/model/covariates/"),
        adata=temp("{base_gex}/annotate/scvi/bcells.h5ad.gz"),
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/scvi.obj_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"],
        ),
        subsample="False",
    conda:
        "scvi-2023"
    shell:
        "python {params.scripts}train_scvi.py {input} -output_file {output.adata} -output_model {output.model} -covariate_genes {params.cell_cycle_genes} -subsample {params.subsample}"


rule merge_vdj:
    input:
        bcell_h5ad="{base_gex}/annotate/scvi/bcells.h5ad.gz",
        vdj_df=ancient(
            "{base_gex}/vdjc/integrated_no_contaminants_with_all_ambient_flags.tsv.gz"
        ),
    output:
        adata="{base_gex}/annotate/vdjc/bcells.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/scvi.obj_preprocess.log",
    run:
        import pandas as pd
        import scanpy as sc

        adata = sc.read_h5ad(input.bcell_h5ad)
        adata.obs_names_make_unique()
        vdj_df = pd.read_table(input.vdj_df)
        print(vdj_df.shape)
        # prepare dfs for merging
        adata.obs["cb"] = (
            adata.obs.reset_index()["index"].str.split("-", expand=True)[0].values
        )
        # sort by umi counts
        vdj_df.sort_values("n_umis", inplace=True, ascending=False)
        # remove duplicates (extra assemblies)
        vdj_df = vdj_df.drop_duplicates(subset=["cb", "sample_uid"], keep="first")
        # select only columns of interest
        vdj_columns = [
            "c_call",
            "v_mismatch",
            "vdj_sequence",
            "v_call",
            "n_umis",
            "n_droplets_vdj",
            "v_identity",
            "lineage_id",
            "cb",
            "sample_uid",
            "locus",
        ]
        df = pd.merge(
            adata.obs,
            vdj_df[vdj_columns],
            how="left",
            suffixes=("", "_y"),
            on=["cb", "sample_uid"],
        )
        adata.obs = df
        adata.obs.set_index("cb", inplace=True)
        adata.obs_names_make_unique()
        adata.write_h5ad(output.adata)


rule subset_bcells:
    input:
        adata="{base_gex}/annotate/vdjc/bcells.h5ad.gz",
    output:
        asc="{base_gex}/outs/asc.h5ad.gz",
        mb="{base_gex}/outs/mb.h5ad.gz",
        nb_other="{base_gex}/outs/nb_other.h5ad.gz",
        all_igh="{base_gex}/outs/all_igh.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/subset_bcells_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/post_cellranger/",
    run:
        import scanpy as sc

        adata = sc.read_h5ad(str(input))
        adata = adata[adata.obs.Immune_All_Low_conf_score > 0.97]
        # add convenience column named celltypist
        adata.obs["celltypist"] = adata.obs["Immune_All_Low_predicted_labels"]
        mb = adata[adata.obs.Immune_All_Low_predicted_labels.str.contains("Memory|Age|Germinal|Proliferative")]
        mb.write_h5ad(output.mb, compression="gzip")
        asc = adata[adata.obs.Immune_All_Low_predicted_labels.str.contains("Plasma")]
        asc.write_h5ad(output.asc, compression="gzip")
        nb_other = adata[
            ~adata.obs.Immune_All_Low_predicted_labels.str.contains("Plasma|Memory|Age|Germinal|Proliferative")
        ]
        nb_other.write_h5ad(output.nb_other, compression="gzip")
        all_igh = adata[~adata.obs.lineage_id.isna()]
        all_igh.write_h5ad(output.all_igh, compression="gzip")
