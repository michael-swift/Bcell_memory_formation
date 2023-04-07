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
parser.add_argument('-covariate_genes', help="tabular data of genes associated with cell cycle, sex, and cell stress")
parser.add_argument('-subsample', default="False", help="whether or not to subsample the input data (for speed of testing different model constructions)")

args = parser.parse_args()
FILENAME = args.input_file
outfile = args.output_file
out_dir = args.output_model
covariate_genes = args.covariate_genes
subsample = args.subsample
print(type(subsample))
covariate_genes_df = pd.read_table(covariate_genes, index_col=0)
nuisance_genes = covariate_genes_df['cc'].to_list() + ["XIST", "FOS", "JUN"]
print("loading data")
adata = sc.read_h5ad(FILENAME)
if subsample == "True":
    sc.pp.subsample(adata, n_obs=50000)
# setup batch column
adata.obs['donor_tissue'] = adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
#setup model
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="donor_tissue")
print(f'training model')
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb", categorical_covariate_keys=["sample_type", "sex", "corr_cycling"])
# train
vae.train(max_epochs=25)
# get latent rep
Z_hat = vae.get_latent_representation()
# calculate umap using latent representation
adata.obsm["X_scVI"] = Z_hat
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, key_added="leiden_scVI", resolution=0.8)
vae.save(out_dir)
add_continuous_covariates = True
if add_continuous_covariates:
    # then copy the expression of each nuisance gene into adata.obs where the key
    # is the gene name
    for g in nuisance_genes:
        adata.obs.loc[:,g] =  sc.get.obs_df(adata, keys = g, layer = 'counts').values
    # finally, remove the nuisance genes from the anndata
    gene_subset = [g for g in adata.var_names if g not in nuisance_genes]
    print(gene_subset)
    adata = adata[:,gene_subset].copy()
    print(adata.shape, "shape of adata after removing nuisance genes")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="donor_tissue", continuous_covariate_keys=nuisance_genes)
    print(f'training model with additional covariates')
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    # train
    vae.train(max_epochs=20)
    # get latent rep
    Z_hat = vae.get_latent_representation()
    # calculate umap using latent representation
    adata.obsm["X_scVI_cov"] = Z_hat
    sc.pp.neighbors(adata, use_rep="X_scVI_cov", n_neighbors=20)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, key_added="leiden_scVI_cov", resolution=0.8)

vae.save(out_dir)
adata.write_h5ad(outfile)
