#!/usr/bin/env python
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from celltypist import models

# single cell processing functions



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

    sc.pl.scatter(adata, x="log10_n_genes", y="log10_total_counts", color="sample_uid")
    if filter_cells == True:
        sc.pp.filter_cells(adata, min_counts=1000)
        sc.pp.filter_cells(adata, min_genes=300)
        sc.pp.filter_cells(adata, max_genes=9000)
        sc.pp.filter_cells(adata, max_counts=89000)
    adata = adata[adata.obs["pct_counts_mt"] < 8]
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
    sc.pl.umap(adata, color="batch")
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
    if batch_correct == True:
        sc.external.pp.bbknn(adata, batch_key="sample_iud")
    else:
        sc.pp.neighbors(adata, n_neighbors=10)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.2)
    sc.pl.umap(adata, color="sample_uid")
    return adata

h5ad, samplesheet = snakemake.input

adata = sc.read_h5ad(h5ad)

adata = perform_qc(adata, filter_cells=True)

adata = add_samplesheet(samplesheet, adata)

print("transforming gene expression")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata, base=2, chunk_size=10000)
sc.pp.highly_variable_genes(adata, n_top_genes=4000)
sc.pl.highly_variable_genes(adata)
adata.raw = adata
sc.pp.scale(adata, max_value=10)

print("annotating with celltypist")
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
sc.tl.rank_genes_groups(adata, groupby="leiden")
sc.tl.rank_genes_groups(adata, groupby="celltypist", key_added="celltypist")

adata.write_h5ad(str(snakemake.output[0]), compression = 'gzip')


# B cells
bcells = adata[adata.obs[anno_1].str.contains("B cells|Plasma")]
adata = None
bcells = recluster(bcells, batch_correct=False)

# write objects

bcells.write_h5ad(str(snakemake.output[1]), compression = 'gzip')
# subsample:
cells = bcells.obs.groupby("sample_uid").sample(frac=0.3).index
bcells = bcells[bcells.obs.index.isin(cells)]

bcells.write_h5ad(str(snakemake.output[2]), compression = 'gzip')
