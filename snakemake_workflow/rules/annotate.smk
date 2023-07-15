rule scvi_all_cells:
    input:
        rules.aggregate_h5ads.output
    output:
        adata="{base_gex}/scvi_all/gex.h5ad.gz",
    log:
        "{base_gex}/logs/annotate/scvi.obj_all.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/gex_preprocess/",
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"],
        ),
        subsample="True",
    conda:
        "scvi-2023"
    shell:
        "python {params.scripts}train_scvi_all.py {input} -output_file {output.adata} -covariate_genes {params.cell_cycle_genes} -subsample {params.subsample}"


rule run_celltypist_by_tissue:
    input:
        rules.scvi_all_cells.output
    output:
        tissue_obj="{base_gex}/annotate/tissue_objs/{tissue}/celltypist.h5ad.gz",
    resources:
        partition="quake,owners",
    params:
    log:
        "{base_gex}/logs/annotate/{tissue}/celltypist_by_tissue.log",
    conda:
        config["workflow_dir"] + "/envs/scanpy.yaml"
    script:
        config["workflow_dir"] + "/scripts/gex_preprocess/run_celltypist.py"

rule annotate_cell_cycle_b_cell:
    """ performs preprocessing and qc, cell cycle assignment, doublet filtering, adds samplesheet info to object"""
    input:
        rules.run_celltypist_by_tissue.output
    output:
        "{base_gex}/annotate/tissue_objs/{tissue}/annotated_processed.h5ad.gz",
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
        config["workflow_dir"] + "/scripts/gex_preprocess/annotate_doublet_cc.py"

rule aggregate_annotated:
    input:
        h5ads=list(set(expand(
            "{base_gex}/annotate/tissue_objs/{tissue}/annotated_processed.h5ad.gz",
            base_gex=base["gex"],
            tissue=tissues,
        ))),
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
        config["workflow_dir"] + "/scripts/gex_preprocess/aggregate_annotated.py"

rule remove_nonb:
    input:
        rules.aggregate_annotated.output.full_h5ad
    output:
        bcells=temp("{base_gex}/annotate/bcells.h5ad.gz"),
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/removenonb_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/gex_preprocess/",
    run:
        import scanpy as sc
        adata = sc.read_h5ad(str(input))
        bcells = adata[(adata.obs.possible_b_cell == True) | (adata.obs.possible_b_cell == "True")]
        bcells.write_h5ad(output.bcells, compression="gzip")

rule dummy_make_vdj:
    "need to connect to vdjc pipeline"
    output:
        "{base_gex}/vdjc/integrated_no_contaminants_with_all_ambient_flags.tsv.gz",
    shell:
        "cat {output}"

rule scvi_bcells:
    input:
        rules.remove_nonb.output
    output:
        model=directory("{base_gex}/annotate/scvi/model/covariates/"),
        adata="{base_gex}/annotate/scvi/bcells.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/scvi.obj_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/gex_preprocess/",
        cell_cycle_genes="{}/resources/cell_cycle/cell_cycle_genes.tab".format(
            config["workflow_dir"],
        ),
        subsample="False",
        layer="counts"
    conda:
        "scvi-2023"
    shell:
        "python {params.scripts}train_scvi.py {input} -output_file {output.adata} -output_model {output.model} -covariate_genes {params.cell_cycle_genes} -subsample {params.subsample} -layer {params.layer}"

rule merge_vdj:
    input:
        bcell_h5ad=rules.scvi_bcells.output.adata,
        vdj_df=ancient(
            "{base_gex}/vdjc/integrated_no_contaminants_with_all_ambient_flags.tsv.gz"
        ),
    output:
        adata=temp("{base_gex}/annotate/vdjc/bcells.h5ad.gz"),
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
            "vdj_sequence",
            "v_call",
            "n_umis",
            "n_droplets_vdj",
            "v_identity",
            "v_mismatch",
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
        asc="{base_gex}/outs/ASC.h5ad.gz",
        mb="{base_gex}/outs/MB.h5ad.gz",
        nb_other="{base_gex}/outs/NB_other.h5ad.gz",
        only_igh="{base_gex}/outs/only_igh.h5ad.gz",
        bcells="{base_gex}/outs/bcells.h5ad.gz",
    resources:
        partition="quake,owners",
    log:
        "{base_gex}/logs/annotate/subset_bcells_preprocess.log",
    params:
        scripts=config["workflow_dir"] + "/scripts/gex_preprocess/",
    run:
        import scanpy as sc
        adata = sc.read_h5ad(str(input))
        # add convenience column named celltypist
        adata.obs["celltypist"] = adata.obs["Immune_All_Low_predicted_labels"]
        MB = adata[
            adata.obs.Immune_All_Low_predicted_labels.str.contains(
                "Memory|Age|Germinal|Proliferative"
            )
        ]
        MB.write_h5ad(output.mb, compression="gzip")
        ASC = adata[adata.obs.Immune_All_Low_predicted_labels.str.contains("Plasma")]
        ASC.write_h5ad(output.asc, compression="gzip")
        NB_other = adata[
            ~adata.obs.Immune_All_Low_predicted_labels.str.contains(
                "Plasma|Memory|Age|Germinal|Proliferative"
            )
        ]
        NB_other.write_h5ad(output.nb_other, compression="gzip")
        only_igh = adata[~adata.obs.lineage_id.isna()]
        only_igh.write_h5ad(output.only_igh, compression="gzip")
        adata.write_h5ad(output.bcells, compression = "gzip")
