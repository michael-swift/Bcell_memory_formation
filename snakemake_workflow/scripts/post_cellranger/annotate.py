import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from celltypist import models
import scrublet as scrub
import gc


###############################
###############################
# single cell preprocessing ###
# functions ###################
###############################
###############################


def perform_qc(adata, filter_cells=True):
    # calculate qc metrics
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True
    )
    sc.pl.highest_expr_genes(adata, save="prefilter")
    # plot qc metrics
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "log1p_total_counts", "total_counts"],
        stripplot=False,
        multi_panel=True,
        save="prefilter",
    )
    if filter_cells == True:
        sc.pp.filter_cells(adata, min_counts=8000)
        sc.pp.filter_cells(adata, min_genes=300)
        sc.pp.filter_cells(adata, max_genes=16000)
        sc.pp.filter_cells(adata, max_counts=300000)
    adata = adata[adata.obs["pct_counts_mt"] < 15]
    # plot results of filtering
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "log1p_total_counts", "total_counts"],
        stripplot=False,
        multi_panel=True,
        save="postfilter",
    )
    sc.pl.highest_expr_genes(adata, save="postfilter")
    return adata


def cluster(adata, batch_correct=False, batch_key="tissue"):
    print("PCA-ing")
    sc.pp.pca(adata)
    print("drawing neighbor graph")
    if batch_correct == True:
        print("batch corrected neighbors")
        sce.pp.bbknn(adata, batch_key=batch_key)
    else:
        sc.pp.neighbors(adata, n_neighbors=20)
    print("UMAP-ing")
    sc.tl.umap(adata)
    print("leiden-ing")
    sc.tl.leiden(adata, resolution=0.2)
    return adata


def add_samplesheet(file, adata):
    samplesheets = pd.read_table(file)
    samplesheets["sample_uid"] = (
        samplesheets["donor"]
        + "_"
        + samplesheets["sample_type"]
        + "_"
        + samplesheets["sample_descriptor"]
    )

    samplesheets.set_index("sample_uid", inplace=True)
    samplesheets = samplesheets[samplesheets.libtype == "gex"]
    for column in samplesheets.columns:
        _dict = samplesheets[column].to_dict()
        adata.obs[column] = adata.obs.sample_uid.map(_dict)
    return adata


def recluster(adata, batch_correct):
    sc.pp.pca(adata)
    print("constructing neighbors graph")
    if batch_correct == True:
        sc.external.pp.bbknn(adata, batch_key="sample_iud")
    else:
        sc.pp.neighbors(adata, n_neighbors=20)
    print("calculating umap")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.2)
    return adata


def setup_figure_outs(tissue):
    output_dir = "{}/figures/QCandAnnotation/gex_{}".format(output_dir, tissue)
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_suffix = ""
    output_formats = [".png", ".svg"]
    sc.settings.figdir = output_dir
    sc.set_figure_params(
        format="pdf",
        transparent=True,
    )

    return output_dir, output_formats, output_suffix


def run_celltypist(adata):
    # create low resolution leiden for majority voting
    sc.tl.leiden(adata, resolution=0.2)
    # High Resolution
    predictions = celltypist.annotate(
        adata, model="Immune_All_Low.pkl", majority_voting=False
    )
    adata = predictions.to_adata(prefix="Immune_All_Low_")
    # not sure why this isn't added automatically
    adata.uns["log1p"] = {"base": np.e}
    # Low Resolution
    predictions = celltypist.annotate(
        adata, model="Immune_All_High.pkl", majority_voting=False
    )
    adata = predictions.to_adata(prefix="Immune_All_High_")
    predictions = celltypist.annotate(
        adata,
        model="Immune_All_Low.pkl",
        majority_voting=True,
        over_clustering=adata.obs.loc[:"leiden"],
    )
    adata = predictions.to_adata(prefix="my_leiden_")
    g = sns.jointplot(
        data=adata.obs,
        x="Immune_All_Low_conf_score",
        y="Immune_All_High_conf_score",
        kind="kde",
    )
    save_figure(g.figure, "confidence_scores_celltypist")
    return adata


def run_scrublet(adata, layers):
    for layer in layers:
        scrub = scr.Scrublet(adata.to_df(layer=layer))
        prediction = "predicted_doublets_{}".format(layer)
        score = "doublet_scores_{}".format(layer)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        adata.obs.loc[:, score] = doublet_scores
        adata.obs.loc[:, prediction] = predicted_doublets
        # modify predicted doublets based on manual score cutoff
        sns.displot(data=adata.obs, x=score)
        plt.yscale("log")
        sns.displot(data=adata.obs, x=score, kind="ecdf")
        adata.obs.loc[:, prediction] = adata.obs[score] > 0.2
        adata.obs.loc[:, prediction] = adata.obs[prediction].astype(str)
        # calculate fraction of leiden cluster that is doublets, conservatively if more than 25 % are doublets, label every cell as doublet associated
        # doing this because scrublet appears to suffer from quite a few false negatives, but setting the cutoff lower is less surgical
        selector = (
            adata.obs.groupby(
                "{}_leiden".format(tissue)
            ).predicted_doublets_umi_counts.value_counts(normalize=True)[:, "True"]
            > 0.25
        )
        exclude_leidens = selector[selector].index
        adata.obs["doublet_associated"] = adata.obs["{}_leiden".format(tissue)].isin(
            exclude_leidens
        )
    # plot layer comparisons
    g = sns.heatmap(
        sc.metrics.confusion_matrix(
            "predicted_doublets_{}".format(layers[0]),
            "predicted_doublets_".format(layers[1]),
            adata.obs,
            normalize=False,
        ),
        annot=True,
    )
    save_figure(g.figure, "doublet_layers_confusion_matrix")
    ## Plot Background Removal Examples
    gene = "IGKC"
    _umi = sc.get.obs_df(adata, keys=gene, layer="umi_counts")
    _bg = sc.get.obs_df(adata, keys=gene, layer="background_removed")
    data = pd.concat([_bg, _umi], axis=1)
    data.columns = ["cbender", "umi"]
    data = data.melt()
    g = sns.displot(data, x="value", hue="variable", kind="ecdf", log_scale=True)
    save_figure(g.figure, "{}_bg_remove".format(gene))

    gene = "AZU1"
    _umi = sc.get.obs_df(adata, keys=gene, layer="umi_counts")
    _bg = sc.get.obs_df(adata, keys=gene, layer="background_removed")
    data = pd.concat([_bg, _umi], axis=1)
    data.columns = ["cbender", "umi"]
    data = data.melt()
    g = sns.displot(data, x="value", hue="variable", kind="ecdf", log_scale=True)
    save_figure(g.figure, "{}_bg_remove".format(gene))
    gene = "IL7R"
    _umi = sc.get.obs_df(adata, keys=gene, layer="umi_counts")
    _bg = sc.get.obs_df(adata, keys=gene, layer="background_removed")
    data = pd.concat([_bg, _umi], axis=1)
    data.columns = ["cbender", "umi"]
    data = data.melt()
    g = sns.displot(data, x="value", hue="variable", kind="ecdf", log_scale=True)
    save_figure(g.figure, "{}_bg_removal_comparison".format(gene))


def assign_technical_categories(adata, tissue):
    """Possible Categories are
    probable_hq_single_b_cell:
    probable_hq_single_not_b_cell:
    possible_b_cell:
    rare_or_bad_q_cell:
    """
    adata.obs["possible_b_cell"] = (
        adata.obs["my_leiden_majority_voting"].str.contains("B cells|Plasma")
        | adata.obs["Immune_All_High_predicted_labels"].str.contains("B cell|Plasma")
        | adata.obs["Immune_All_Low_predicted_labels"].str.contains("B cell|Plasma")
    )
    # change categorical to boolian
    bool_mapper = {"True": True, "False": False}
    adata.obs["predicted_doublets_umi_counts"] = (
        adata.obs["predicted_doublets_umi_counts"].astype(str).map(bool_mapper)
    )
    adata.obs["probable_hq_single_b_cell"] = (
        adata.obs["possible_b_cell"]
        & ~adata.obs["predicted_doublets_umi_counts"]
        & (~adata.obs["{}_leiden".format(tissue)].astype(int) <= 16)
        & (~adata.obs["doublet_associated"])
        & adata.obs["Immune_All_High_predicted_labels"].str.contains("B cell|Plasma")
        & adata.obs["Immune_All_Low_predicted_labels"].str.contains("B cell|Plasma")
    )

    adata.obs["probable_hq_single_not_b_cell"] = (
        ~adata.obs["possible_b_cell"]
        & ~adata.obs["predicted_doublets_umi_counts"]
        & (~adata.obs["{}_leiden".format(tissue)].astype(int) <= 16)
        & (~adata.obs["doublet_associated"])
    )

    adata.obs["rare_or_bad_q_cell"] = (
        adata.obs["{}_leiden".format(tissue)].astype(int) > 16
    )


def annotate_cell_cyle(adata, cell_cycle_genes):
    # develop a crude way to label cycling and non-cycling then test features by t-test
    genes_contributing_to_score = 30
    sc.tl.score_genes(
        adata,
        gene_list=cell_cycle_genes.loc[:genes_contributing_to_score, "cc"],
        score_name="correlation_cycling",
    )
    sc.tl.score_genes(
        adata,
        gene_list=cell_cycle_genes.loc[:genes_contributing_to_score, "anti_cc"],
        score_name="anticorrelation_cycling",
    )
    adata.obs["corr_cycling"] = adata.obs["correlation_cycling"] > 0.3
    adata.obs["anticorr_cycling"] = adata.obs["anticorrelation_cycling"] > 0.4
    cycling_mapper = {True: "True", False: "False"}
    use_non_cycling = False
    if use_non_cycling:
        adata.obs["cycling"] = adata.obs["anticorr_cycling"].map(cycling_mapper)
    else:
        adata.obs["cycling"] = adata.obs["corr_cycling"].map(cycling_mapper)
    fig, axs = plt.subplots(1, 2, sharey=True)
    fig.tight_layout(pad=2.0)
    x, y = adata.obs["anticorrelation_cycling"], adata.obs["correlation_cycling"]
    sns.histplot(x, ax=axs[0])
    sns.histplot(y, ax=axs[1])
    plt.vlines(x=0.1, ymax=10000, ymin=0, linestyles="dotted", color="r")
    plt.yscale("log")
    save_figure(fig, "{}_cell_cycle_scoring".format(tissue))
    fig, axs = plt.subplots(1, 2, sharey=True)
    fig.tight_layout(pad=2.0)
    x, y = adata.obs["anticorrelation_cycling"], adata.obs["correlation_cycling"]
    sns.histplot(x, ax=axs[0])
    sns.histplot(y, ax=axs[1])
    plt.vlines(x=0.1, ymax=10000, ymin=0, linestyles="dotted", color="r")
    plt.yscale("log")
    save_figure(fig, "{}_cell_cycle_scoring".format(tissue))
    return adata


def run_gex_label_routine(adata, tissue):
    output_dir, output_formats, output_suffix = setup_figure_outs(tissue=tissue)

    def save_figure(
        fig,
        name,
        output_dir=output_dir,
        output_suffix=output_suffix,
        output_formats=output_formats,
        savefig_args=savefig_args,
    ):
        for output_format in output_formats:
            fig.savefig(
                output_dir + "/" + name + output_suffix + output_format, **savefig_args
            )

    # remove IGH and IGL genes
    # and genes associated with dissociation for clustering analysis
    # just from variable gene definition
    adata_tissue.var.loc[
        adata_tissue.var.index.str.contains("IGH|IGL|IGK|FOS|JUN|HSP|RPL"),
        "highly_variable",
    ] = False
    # ad hoc exclusion of weird samples, this sample is so overloaded I don't think we can analyze the GEX data
    adata_tissue = adata_tissue[adata_tissue.obs.sample_uid != "TBd3_fresh_B200"]
    adata_tissue = cluster(adata_tissue, batch_correct=False)
    filter_low_abundance_cell_groups = False
    cell_group = "predicted_labels"
    if filter_low_abundance_cell_groups:
        select = adata_tissue.obs[cell_group].value_counts() > (
            adata_tissue.obs.shape[0] / 1000
        )
        adata_tissue = adata_tissue[
            adata_tissue.obs[cell_group].isin(select[select == True].index)
        ]
    #### Celltypist ####
    adata_tissue = run_celltypist(adata_tissue)
    #### Scrublet ####
    # doublet removal and scoring
    adata_tissue = run_scrublet(adata_tissue, layers=["counts", "cellbender_counts"])
    # Algo for definitely being a B cell
    assign_technical_categories(adata_tissue, tissue)
    # Cell cycle Annotation
    adata_tissue = annotate_cell_cycle(adata_tissue, cell_cycle_genes)

    return adata_tissue


###############################
#### execution ################
###############################

h5ad, samplesheet = snakemake.input
cell_cycle_genes = snakemake.params.cell_cycle_genes
tissue_dir = snakemake.output.tissue_dir
tissues = snakemake.wildcards.tissues
adata = sc.read_h5ad(h5ad)
adata.var_names_make_unique()
adata.obs_names_make_unique()
subsample = True
if subsample:
    adata = adata[adata.obs.index.isin(adata.obs.sample(n=50000, replace=False).index)]
    print("subsampling data to ", adata.shape[0], " cells")
# QC and metadata integration
adata = perform_qc(adata, filter_cells=True)
adata = add_samplesheet(samplesheet, adata)
# data transform
print("transforming gene expression")
print("normalize per 10K counts")
sc.pp.normalize_total(adata, target_sum=1e4)
print("log base 2")
sc.pp.log1p(adata, chunk_size=10000, base=2)
print("top highly variable genes")
adata.raw = adata
# perform gex labeling routine
tissue = snakemake.wildcards.tissue
print(tissue)
gex_label_to_write = [
    "sample_uid",
    "Immune_All_Low_predicted_labels",
    "Immune_All_High_predicted_labels",
    "Immune_All_High_conf_score",
    "Immune_All_Low_conf_score",
    "doublet_scores_umi_counts",
    "predicted_doublets_umi_counts",
    "n_genes",
    "total_counts",
    "total_counts_mt",
    "{}_leiden".format(tissue),
    "tissue",
    "probable_hq_single_b_cell",
    "possible_b_cell",
    "probable_hq_single_not_b_cell",
    "rare_or_bad_q_cell",
]
adata = adata[adata.obs["tissue"] == tissue]
print(adata.shape)
adata_tissue = run_gex_label_routine(adata, tissue)
# write obs matrix for vdj integration
adata_tissue.write_h5ad(
    "{}/{}_annotated_processed.h5ad.gz".format(tissue_dir, tissue),
    compression="gzip",
)
df = adata_tissue.obs[gex_labels_to_write]
df.to_csv("{}/{}_gex_labels.tab".format(tissue_dir, tissue), sep="\t")
# perform garbage collection:
del adata_tissue
gc.collect()
