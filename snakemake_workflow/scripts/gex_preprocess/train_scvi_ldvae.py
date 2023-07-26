#!/usr/bin/env python
import numpy as np
import pandas as pd
import scvi
import argparse
import scanpy as sc
import torch

# I was told this speeds up the calculations
torch.set_float32_matmul_precision("medium")
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
parser.add_argument(
    "-layer",
    default="counts",
    help="which data layer to use e.g. 'counts' or 'cellbender_counts' or not to subsample the input data (for speed of testing different model constructions)",
)

parser.add_argument(
    "-latent_rep",
    default="scvi_latent.csv",
    help="the latent representation calculated by scvi, useful for adding to other anndatas without needing to load the whole h5ad",
)

parser.add_argument(
    "-loadings",
    default="ldvae_loadings.csv",
    help="the latent representation calculated by scvi, useful for adding to other anndatas without needing to load the whole h5ad",
)

args = parser.parse_args()

FILENAME = args.input_file
outfile = args.output_file
out_dir = args.output_model
subsample = args.subsample
layer = args.layer
loadings = args.loadings
latent_rep = args.latent_rep
print(type(subsample))
print("loading data")
adata = sc.read_h5ad(FILENAME)
# move to params at somepoint
batch_key = 'tissue'
add_linear_decoded = True
if subsample == 'True':
   sc.pp.subsample(adata, n_obs = 10000)
if add_linear_decoded:
    # setup linear
    layer = "counts"
    scvi.model.LinearSCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
    model = scvi.model.LinearSCVI(adata, n_latent=8)
    # train
    model.train(max_epochs=150, plan_kwargs={"lr":5e-3}, check_val_every_n_epoch=10)
    # get latent rep
    Z_hat = model.get_latent_representation()
    # add latent rep to obsm
    adata.obsm["X_scVI_bcells_LDVAE"] = Z_hat
    # extract loadings and also add as a varm:
    loadings = model.get_loadings()
    print(loadings.head())
    adata.varm["LDVAE_loadings"] = loadings
print("writing output")
adata.write_h5ad(outfile)
adata.obsm['X_scVI_bcells_LDVAE'].to_csv(latent_rep)
adata.varm["LDVAE_loadings"].to_csv(loadings)
