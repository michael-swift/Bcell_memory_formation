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

args = parser.parse_args()

FILENAME = args.input_file
outfile = args.output_file
out_dir = args.output_model
subsample = args.subsample
layer = args.layer
print(type(subsample))
print("loading data")
adata = sc.read_h5ad(FILENAME)
if subsample == "True":
    sc.pp.subsample(adata, n_obs=100000)
# setup batch column
adata.obs["donor_tissue"] = (
    adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
)
print(adata.obs['donor_tissue'].value_counts())

batch_key = 'tissue'


# setup scVI model:
scvi.model.JaxSCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
vae = scvi.model.JaxSCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# train
vae.train(max_epochs=150)
# get latent rep
Z_hat = vae.get_latent_representation()
# calculate umap using latent representation
adata.obsm["X_scVI_b_cells"] = Z_hat
add_continuous = True
if add_continuous:
    adata_cond = adata.copy()
    # cell cycle genes, IGH genes, and dissociation gene 
    # not exhaustive but bc they covary strongly these three genes appear to take care of the cell cyle effect 
    cell_cycle_genes = ["MKI67", "RRM2", "TK1"]
    nuisance_genes = ["FOS", "XIST", "JUN", "DNAJB1", "HSPA1B"] 
    ighl_genes = ["IGHD", "IGHA1", "IGHA2", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHE", "IGLC3", "IGLC1", "IGLC2","IGKC"]
    # could consider removing IGHV genes here as well but they seem to help identify ambient clusters
    for g in nuisance_genes + ighl_genes + cell_cycle_genes:
        exp = adata_cond[:, g].X.toarray()
        adata_cond.obs[g] = exp.copy()
    # finally, remove the nuisance genes from the anndata
    gene_subset = [g for g in adata.var_names if g not in nuisance_genes]
    adata_cond = adata_cond[:, gene_subset].copy()

    scvi.model.SCVI.setup_anndata(
        adata_cond,
        layer=layer,
        batch_key=batch_key,
        continuous_covariate_keys=nuisance_genes,
    )
    vae = scvi.model.SCVI(adata_cond, n_layers=2, n_latent=20, gene_likelihood="nb")
    # train
    vae.train(max_epochs=150)
    # get latent rep
    Z_hat = vae.get_latent_representation()
    # calculate umap using latent representation
    adata.obsm["X_scVI_bcells_nuisance"] = Z_hat
    vae.save(out_dir)

add_linear_decoded = True
if add_linear_decoded:
    # setup linear
    layer = "counts"
    scvi.model.LinearSCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
    model = scvi.model.LinearSCVI(adata, n_latent=20)
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

adata.write_h5ad(outfile)