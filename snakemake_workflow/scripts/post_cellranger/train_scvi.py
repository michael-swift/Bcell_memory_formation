#!/usr/bin/env python
import numpy as np
import pandas as pd
import scvi
import argparse
import scanpy as sc

# read h5ad
parser = argparse.ArgumentParser(description="run scvi")
parser.add_argument("input_file", help="h5ad anndata")
parser.add_argument("-output_file", help="name of output file")
parser.add_argument("-covariate_genes", help=".tsv with genes, not currently used")
parser.add_argument(
    "-output_model",
    default=".",
    help="model output directory (default: working directory)",
)
parser.add_argument(
    "-subsample",
    default="False",
    help="whether or not to subsample the input data (for speed of testing different model constructions)",
)

args = parser.parse_args()
FILENAME = args.input_file
outfile = args.output_file
out_dir = args.output_model
subsample = args.subsample
print(type(subsample))
print("loading data")
adata = sc.read_h5ad(FILENAME)

layer = "counts"
if subsample == "True":
    sc.pp.subsample(adata, n_obs=50000)
# setup batch column
adata.obs["donor_tissue"] = (
    adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
)

# setup scVI model:
scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key="donor_tissue")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# train
vae.train(max_epochs=150)
# get latent rep
Z_hat = vae.get_latent_representation()
# calculate umap using latent representation
adata.obsm["X_scVI_raw"] = Z_hat
add_continuous = True
if add_continuous:
    adata_cond = adata.copy()
    nuisance_genes = ["FOS", "JUN", "DNAJB1", "MKI67", "RRM2", "TK1", "HSPA1B"]
    # then copy the expression of each nuisance gene into adata.obs where the key
    # is the gene name
    for g in nuisance_genes:
        exp = adata_cond[:, g].X.toarray()
        adata_cond.obs[g] = exp.copy()
    # finally, remove the nuisance genes from the anndata
    gene_subset = [g for g in adata.var_names if g not in nuisance_genes]
    adata_cond = adata_cond[:, gene_subset].copy()

    scvi.model.SCVI.setup_anndata(
        adata_cond,
        layer=layer,
        batch_key="donor_tissue",
        continuous_covariate_keys=nuisance_genes,
    )
    vae = scvi.model.SCVI(adata_cond, n_layers=2, n_latent=30, gene_likelihood="nb")
    # train
    vae.train(max_epochs=150)
    # get latent rep
    Z_hat = vae.get_latent_representation()
    # calculate umap using latent representation
    adata.obsm["X_scVI_raw_cont"] = Z_hat
    vae.save(out_dir)
add_cellbender_model = False
if add_cellbender_model:
    # setup vanilla model:
    layer = "cellbender_counts"
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key="donor_tissue")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    # train
    vae.train(max_epochs=150)
    # get latent rep
    Z_hat = vae.get_latent_representation()
    # calculate umap using latent representation
    adata.obsm["X_scVI_cellbender"] = Z_hat

adata.write_h5ad(outfile)
