import scanpy as sc
import pandas as pd

file = str(snakemake.input.cr)
print(file)
sample_uid = snakemake.wildcards.sample_uid
min_genes_per_CB = snakemake.params.min_genes
min_counts_per_CB = snakemake.params.min_counts

# read 10X Genomics raw file
adata = sc.read_10x_h5(file)
adata.obs['sample_uid'] = sample_uid

sc.pp.filter_cells(adata, min_genes = min_genes_per_CB)
sc.pp.filter_cells(adata, min_counts = min_counts_per_CB)
# make index names unique
#adata.obs_names_make_unique()
#adata.var_names_make_unique()

adata.write_h5ad(str(snakemake.output))
print("Done")




