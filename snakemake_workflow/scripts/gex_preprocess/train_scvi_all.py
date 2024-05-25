#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Run scVI")
    parser.add_argument("input_file", help="h5ad anndata")
    parser.add_argument("-output_file", help="Name of output file")
    parser.add_argument(
        "-output_model",
        default=".",
        help="Model output directory (default: working directory)"
    )
    parser.add_argument(
        "-covariate_genes",
        help="Tabular data of genes associated with cell cycle, sex, and cell stress"
    )
    parser.add_argument(
        "-subsample",
        default="False",
        help="Whether or not to subsample the input data (for speed of testing different model constructions)"
    )

    args = parser.parse_args()
    filename = args.input_file
    outfile = args.output_file
    out_dir = args.output_model
    covariate_genes = args.covariate_genes
    subsample = args.subsample

    print(type(subsample))
    print("Loading data")
    adata = sc.read_h5ad(filename)
    layer = "counts"

    if subsample == "True":
        sc.pp.subsample(adata, n_obs=100000)

    # Setup batch column
    adata.obs["donor_tissue"] = (
        adata.obs["donor"].astype(str) + "_" + adata.obs["tissue"].astype(str)
    )

    # Setup vanilla model
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key="donor_tissue")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

    # Train
    vae.train(max_epochs=150)

    # Get latent representation
    z_hat = vae.get_latent_representation()

    # Calculate UMAP using latent representation
    adata.obsm["X_scVI_all"] = z_hat

    add_cellbender_model = False
    if add_cellbender_model:
        # Setup vanilla model
        layer = "cellbender_counts"
        scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key="donor_tissue")
        vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

        # Train
        vae.train(max_epochs=150)

        # Get latent representation
        z_hat = vae.get_latent_representation()

        # Calculate UMAP using latent representation
        adata.obsm["X_scVI_cellbender_all"] = z_hat

    adata.write_h5ad(outfile)

if __name__ == "__main__":
    main()
