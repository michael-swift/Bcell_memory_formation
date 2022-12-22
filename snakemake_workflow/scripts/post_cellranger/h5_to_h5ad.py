import scanpy as sc

tenX = str(snakemake.input.cr)
cellbender = str(snakemake.input.cb)
print(file_name)
sample_uid = snakemake.wildcards.sample_uid
print(sample_uid)
min_genes_per_CB = snakemake.params.min_genes
min_counts_per_CB = snakemake.params.min_counts

# 10X Genomics File:
adata = sc.read_10x_h5(tenX)
adata.obs['sample_uid'] = sample_uid
sc.pp.filter_cells(adata, min_genes = min_genes_per_CB)
sc.pp.filter_cells(adata, min_counts = min_counts_per_CB)
adata.write_h5ad(str(snakemake.output[0]))

# CellBender File
adata = sc.read_10x_h5(cellbender)
adata.obs['sample_uid'] = sample_uid
sc.pp.filter_cells(adata, min_genes = min_genes_per_CB)
sc.pp.filter_cells(adata, min_counts = min_counts_per_CB)
adata.write_h5ad(str(snakemake.output[1]))
print("Done")
