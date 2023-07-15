import scanpy as sc
import celltypist
import numpy as np
import pandas as pd
import scanpy.external as sce
from celltypist import models
import pathlib
import seaborn as sns
import anndata as ad
import gc


def perform_qc(adata, filter_cells=True, min_genes=300, min_counts=1000):
    # calculate qc metrics
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True
    )
    if filter_cells == True:
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.filter_cells(adata, max_genes=20000)
    adata = adata[adata.obs["pct_counts_mt"] < 13.5]
    return adata

def run_celltypist(adata):
    """
    This function annotates a single-cell dataset using CellTypist models.

    Parameters
    ----------
    adata : anndata.AnnData
        The annotated data matrix of shape n_obs x n_vars. Rows correspond
        to cells and columns to genes.

    Returns
    -------
    adata : anndata.AnnData
        The input adata object, updated with new annotations from the
        CellTypist models. Annotations include cell type labels and
        corresponding confidence scores.
    """

    # High resolution annotation using "Immune_All_Low.pkl" model
    # Setting majority_voting to False to avoid over-clustering
    predictions = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=False)

    # Update adata with the new annotations
    adata = predictions.to_adata(prefix="Immune_All_Low_")

    # Explicitly add log1p parameter to adata.uns
    # This step might be necessary for downstream analysis in scanpy
    adata.uns["log1p"] = {"base": np.e}
    # Low resolution annotation using "Immune_All_High.pkl" model
    # Again, set majority_voting to False to avoid over-clustering
    predictions = celltypist.annotate(adata, model="Immune_All_High.pkl", majority_voting=False)
    # Update adata with the new annotations
    adata = predictions.to_adata(prefix="Immune_All_High_")
    # Additional annotation using "Immune_All_Low.pkl" model
    # Here, majority_voting is set to True and over_clustering is enabled
    predictions = celltypist.annotate(
        adata,
        model="Immune_All_Low.pkl",
        majority_voting=True,
    )
    # Update adata with the new annotations
    adata = predictions.to_adata(prefix="majority_voting_low_")
    return adata

def perform_preprocess(adata):
    print("transforming gene expression")
    print("normalize per 10K counts")
    sc.pp.normalize_total(adata, target_sum=1e4)
    print("log base e")
    sc.pp.log1p(adata, chunk_size=10000)
    # remove this sample (over loaded so much it GEX is uninterpretable
    adata = adata[adata.obs.sample_uid != "TBd3_fresh_B200"]
    print("top highly variable genes")
    sc.pp.highly_variable_genes(adata, n_top_genes=1000)
    # remove IGH and Dissociation associated genes from highly variable
    adata.var.loc[adata.var.index.str.contains("IGH|IGL|IGK|FOS|JUN|HSP|RPL"), "highly_variable"
] = False
    adata.raw = adata
    return adata

# load data
h5ad = str(snakemake.input)
output_file = str(snakemake.output)
adata_full = sc.read_h5ad(h5ad, backed = 'r')
tissue = snakemake.wildcards.tissue
print(tissue)
print("shape of adata", adata_full.shape)
adata = adata_full[adata_full.obs["tissue"] == tissue]
adata = adata.to_memory()
del adata_full
gc.collect()
print("shape of adata post filter", adata.shape)
adata.X = adata.layers["counts"]
adata.var_names_make_unique()
adata.obs_names_make_unique()
print(adata.obs.groupby("donor").tissue.value_counts())


adata = perform_qc(adata, filter_cells=True, min_genes=300, min_counts=1000)
adata = perform_preprocess(adata)
adata = run_celltypist(adata)
adata.write_h5ad(output_file)
