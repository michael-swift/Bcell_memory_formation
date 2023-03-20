#!/usr/bin/env python
import numpy as np
import pandas as pd
import scvi
import argparse
import scanpy as sc
# read h5ad
parser = argparse.ArgumentParser(description='run scvi')
parser.add_argument('input_file', help="h5ad anndata")
parser.add_argument('-output_file', help="name of output file")
parser.add_argument('-output_model', default='.', help="model output directory (default: working directory)")


args = parser.parse_args()
FILENAME = args.input_file
outfile = args.output_file
out_dir = args.output_model

print("loading data")
adata = sc.read_h5ad(FILENAME)
# setup batch column
adata.obs['donor_tissue'] = adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
#setup model

print(f'setting up model with batch key')
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="donor_tissue")
print(f'training model')
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
# train
vae.train(max_epochs=150)
# get latent rep
Z_hat = vae.get_latent_representation()
# calculate umap using latent representation
adata.obsm["X_scVI"] = Z_hat
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.8)

vae.save(out_dir)
adata.write_h5ad(outfile)
