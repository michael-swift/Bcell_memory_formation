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
parser.add_argument(
    "-output_model",
    default=".",
    help="model output directory (default: working directory)",
)
parser.add_argument(
    "-covariate_genes",
    help="tabular data of genes associated with cell cycle, sex, and cell stress",
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
covariate_genes = args.covariate_genes
subsample = args.subsample
print(type(subsample))
print("loading data")
adata = sc.read_h5ad(FILENAME)
layer = "counts"
if subsample == "True":
    sc.pp.subsample(adata, n_obs=100000)

# setup batch column
adata.obs["donor_tissue"] = (
    adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
)

# setup vanilla model:
scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key="donor_tissue")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# train
vae.train(max_epochs=150)
# get latent rep
Z_hat = vae.get_latent_representation()
# calculate umap using latent representation
adata.obsm["X_scVI_all"] = Z_hat

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
    adata.obsm["X_scVI_cellbender_all"] = Z_hat
adata.write_h5ad(outfile)
