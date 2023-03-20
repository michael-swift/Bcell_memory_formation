import scanpy as sc
import anndata as ad
import numpy as np

# get files from snake make
files = snakemake.input

# params from snakemake

min_genes = int(snakemake.params.min_genes)
min_counts = int(snakemake.params.min_counts)
filter_cells = snakemake.params.filter_cells

def load_and_filter(filename, min_genes=min_genes, min_counts=min_counts):
    if "h5ad" in filename:
        adata = sc.read_h5ad(filename)
    else:
        adata = sc.read_10x_h5(filename)
    if filter_cells == True:
        sc.pp.filter_cells(adata, min_genes = min_genes)
        sc.pp.filter_cells(adata, min_counts = min_counts)
    adata.obs_names_make_unique()
    adata.var_names_make_unique(join='-')
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
