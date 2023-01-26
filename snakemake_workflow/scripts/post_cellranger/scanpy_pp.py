#!/usr/bin/env python
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from celltypist import models

###############################
###############################
# single cell preprocessing ##
# functions
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
        sc.pp.filter_cells(adata, min_counts=1500)
        sc.pp.filter_cells(adata, min_genes=500)
        sc.pp.filter_cells(adata, max_genes=12000)
        sc.pp.filter_cells(adata, max_counts=150000)
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
    samplesheets = samplesheets[samplesheets.libtype == 'gex']
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

###############################
#### execution ################
###############################

h5ad, samplesheet = snakemake.input

adata = sc.read_h5ad(h5ad)
adata.var_names_make_unique()
adata.obs_names_make_unique()

subsample = False
if subsample:
    adata = adata[adata.obs.index.isin(adata.obs.sample(n=20000, replace = False).index)]
    print(adata.shape)

adata = perform_qc(adata, filter_cells=True)
adata = add_samplesheet(samplesheet, adata)
print("transforming gene expression")
adata.layers['umi_counts'] = adata.X.copy()
print("normalize per 10K counts")
sc.pp.normalize_total(adata, target_sum=1e4)
print("log base 2")
sc.pp.log1p(adata, chunk_size=10000, base = 2)
print("top highly variable genes")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key = "tissue")

# remove IGH and IGL variable genes from highly variable genes for clustering analysis 
adata.var.loc[adata.var.index.str.contains("IGHV|IGLV|IGKV"), 'highly_variable'] = False

# set raw as the normalized log counts
adata.raw = adata

# cluster data before celltypist
adata = recluster(adata, batch_correct = False)
print("annotating with celltypist")

# celltypist annotation
anno = "celltypist"

# Download all the available models.
models.download_models()
# Provide the input as an `AnnData`.
# majority voting forces construction of a neighborhood graph
majority_voting = True
if majority_voting:
    predictions = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
    print('calculating UMAP')
    sc.tl.umap(adata)
else:
    predictions = celltypist.annotate(adata, model = "Immune_All_Low.pkl")

# Get an `AnnData` with predicted labels embedded into the cell metadata columns.
adata = predictions.to_adata()
# this line allows to write h5ad, otherwise h5py throws errors
adata.obs = adata.obs.astype(str, errors = "ignore")
# write h5ads
out_prefix = str(snakemake.output[0]).split('.')[0]
for tissue in adata.obs.tissue.unique():
    print('writing per tissue {} h5ads'.format(tissue))
    sub_adata = adata[adata.obs.tissue == tissue]
    sub_adata.write_h5ad("{}_{}.h5ad.gz".format(out_prefix, str(tissue)), compression = "gzip")

print('writing full-h5ad')
adata.write_h5ad(str(snakemake.output[0]), compression = "gzip")

sc.pp.subsample(adata, fraction=0.2, copy=False)

adata.write_h5ad("{}_{}.h5ad.gz".format(out_prefix, "subsampled"), compression = "gzip")
