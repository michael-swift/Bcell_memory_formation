# import scanpy
import scanpy as sc

# load the data
adata = sc.read_10x_h5('tiny_10x_pbmc_filtered.h5', genome='background_removed')

import tables
import numpy as np

z = []
with tables.open_file('tiny_10x_pbmc_filtered.h5') as f:
    print(f)  # display the structure of the h5 file
    z = f.root.background_removed.latent_gene_encoding.read()  # read latents

np.savetxt('tiny_10x_pbmc_latent_gene_expression.csv', z, delimiter=',')

# load the latent representation into a new slot called 'X_cellbender'
adata.obsm['X_cellbender'] = z

# perform louvain clustering using the cellbender latents and cosine distance
sc.pp.neighbors(adata, use_rep='X_cellbender', metric='cosine')
sc.pp.louvain(adata)
