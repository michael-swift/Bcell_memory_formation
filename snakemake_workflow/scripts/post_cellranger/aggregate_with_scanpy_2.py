import scanpy as sc
import anndata as ad
import tables
import numpy as np

# get files from snake make
cr_files = snakemake.input.cr
cb_files = snakemake.input.cb
cb_cell_calls = snakemake.input.cb_calls

print(cr_files)
print(cb_files)
print(cb_cell_calls)

"""
# params from snakemake
min_genes = int(snakemake.params.min_genes)
min_counts = int(snakemake.params.min_counts)
filter_cells = snakemake.params.filter_cells

def load_and_filter(filename, min_genes=min_genes, min_counts=min_counts):
    if "h5ad" in filename:
        adata = sc.read_h5ad(filename)
    else:
        adata = sc.read_10x_h5(filename)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    # flexibly parse this later for sample_uid
    adata.obs['sample_uid'] = str(filename)
    if filter_cells == True:
        sc.pp.filter_cells(adata, min_genes = min_genes)
        sc.pp.filter_cells(adata, min_counts = min_counts)
    return adata

adata_list = []
# construct aggregated object
for i, file_name in enumerate(files):
    adata_list.append(load_and_filter(file_name))
    print(i, " anndatas concatenated")

adata = ad.concat(adata_list)
print(str(snakemake.output))
adata.write_h5ad(str(snakemake.output))
print("Done!!")



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
"""
