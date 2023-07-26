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

args = parser.parse_args()

FILENAME = args.input_file
outfile = args.output_file
out_dir = args.output_model
subsample = args.subsample
layer = args.layer
latent_rep = args.latent_rep
print(type(subsample))
print("loading data")
adata = sc.read_h5ad(FILENAME)

if subsample == "True":
    sc.pp.subsample(adata, n_obs=10000)
# setup batch column
adata.obs["donor_tissue"] = (
    adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
)
print(adata.obs['donor_tissue'].value_counts())

batch_key = 'tissue'

add_continuous = True
if add_continuous:
    # cell cycle genes, IGH genes, and dissociation gene 
    # not exhaustive but bc they covary strongly these three genes appear to take care of the cell cyle effect 
    cell_cycle_genes = ["MKI67", "RRM2", "TK1"]
    nuisance_genes = ["FOS", "XIST", "JUN", "DNAJB1", "HSPA1B"] + cell_cycle_genes
    # could consider removing IGHV genes here as well but they seem to help identify ambient clusters
    for g in nuisance_genes:
        exp = adata[:, g].X.toarray()
        adata.obs[g] = exp.copy()
    # finally, remove the nuisance genes from the anndata
    gene_subset = [g for g in adata.var_names if g not in nuisance_genes]
    adata = adata[:, gene_subset].copy()

    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        continuous_covariate_keys=nuisance_genes,
    )
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=10, gene_likelihood="nb")
    # train
    vae.train(max_epochs=150)
    # get latent rep
    Z_hat = vae.get_latent_representation()
    # calculate umap using latent representation
    adata.obsm["X_scVI_bcells_nuisance"] = Z_hat
    vae.save(out_dir)
adata.write_h5ad(outfile)

pd.DataFrame(adata.obsm['X_scVI_bcells_nuisance']).to_csv(latent_rep)
