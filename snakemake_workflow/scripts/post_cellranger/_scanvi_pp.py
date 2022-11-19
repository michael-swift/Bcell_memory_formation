#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from celltypist import models
import scvi

# single cell preprocessing functions
def perform_qc(adata, filter_cells=False):
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
    adata.obs["log10_total_counts"] = np.log10(adata.obs["total_counts"] + 1)
    adata.obs["log10_n_genes"] = np.log10(adata.obs["n_genes"] + 1)
    if filter_cells == True:
        sc.pp.filter_cells(adata, min_counts=100)
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_cells(adata, max_genes=10000)
        sc.pp.filter_cells(adata, max_counts=100000)
    adata = adata[adata.obs["pct_counts_mt"] < 12]
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
        sc.pp.neighbors(adata, n_neighbors=10)
    print("calculating umap")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.2)
    return adata

h5ad, samplesheet = snakemake.input

adata = sc.read_h5ad(h5ad)

adata = perform_qc(adata, filter_cells=True)

adata = add_samplesheet(samplesheet, adata)

adata.var_names_make_unique()
adata = adata[adata.obs.sample(frac=0.2, replace = False).index]
print("transforming gene expression")
adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="sample_uid"
)


print("setting up object")
scvi.data.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["tissue", "donor"],
    continuous_covariate_keys=["percent_mito"]
)

print("training the model SCVI")
model = scvi.model.SCVI(adata)

model.train()


latent = model.get_latent_representation()

adata.obsm["X_scVI"] = latent

denoised = model.get_normalized_expression(adata, library_size=1e4)

print(denoised.iloc[:5, :5], "denoised data")

adata.layers["scvi_normalized"] = model.get_normalized_expression(
    library_size=10e4
)

# celltypist annotation
anno_1 = "celltypist"
# Download all the available models.
models.download_models()
# Provide the input as an `AnnData`.
predictions = celltypist.annotate(
    adata, model="Immune_All_Low.pkl", majority_voting=True
)
# Celltypist Annotations to bcells object
adata.obs[anno_1] = predictions.predicted_labels.majority_voting

print("clustering and umaping")
adata = cluster(adata, batch_correct=False)
adata.write_h5ad(str(snakemake.output[0]), compression = 'gzip')

# B cells
bcells = adata[adata.obs[anno_1].str.contains("B cells|Plasma")]
adata = None
bcells = recluster(bcells, batch_correct=False)

# write objects

bcells.write_h5ad(str(snakemake.output[1]), compression = 'gzip')

