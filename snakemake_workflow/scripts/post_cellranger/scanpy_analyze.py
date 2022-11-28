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
        sc.pp.filter_cells(adata, min_counts=1000)
        sc.pp.filter_cells(adata, min_genes=300)
        sc.pp.filter_cells(adata, max_genes=9000)
        sc.pp.filter_cells(adata, max_counts=89000)
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

###############################
#### execution ################
###############################

h5ad, samplesheet = snakemake.input

adata = sc.read_h5ad(h5ad)
adata.obs['sample_uid'] = adata.obs['sample_uid'].str.rsplit('/', expand = True).iloc[:,10]
adata = perform_qc(adata, filter_cells=True)
adata = add_samplesheet(samplesheet, adata)
adata.var_names_make_unique()
adata.obs_names_make_unique()

print("transforming gene expression")
adata.layers['raw_counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata, base=2, chunk_size=10000)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pl.highly_variable_genes(adata)

# set raw as the normalized log counts
adata.raw = adata
print("annotating with celltypist")
# celltypist annotation
anno = "celltypist"
# Download all the available models.
models.download_models()
# Provide the input as an `AnnData`.
predictions = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
# Celltypist Annotations to bcells object
adata.obs[anno] = predictions.predicted_labels.majority_voting
print(adata.obs[anno])
print('writing h5ad')
adata.write_h5ad(str(snakemake.output[0]))
print("writing h5ads by sample_uid")

for uid in adata.obs.sample_uid.unique():
    sub_adata = adata[adata.obs.sample_uid == uid]
    sub_adata.write_h5ad("{}/{}.h5ad".format(str(snakemake.params.outdir), uid))

write_flat_files = False

if write_flat_files:
    # write flat files since read AnnData isn't working for scTK
    t = adata.layers['raw_counts'].toarray()
    pd.DataFrame(data = t, index = adata.obs_names, columns = adata.var_names).to_csv(str(snakemake.output[1]), sep = '\t')
    adata.obs.to_csv(str(snakemake.output[2]), sep = '\t')
    adata.var.to_csv(str(snakemake.output[3]), sep = '\t')
