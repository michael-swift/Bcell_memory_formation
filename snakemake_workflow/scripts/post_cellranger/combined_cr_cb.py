import scanpy as sc
import pandas as pd

cellranger = str(snakemake.input.cr)
cellbender = str(snakemake.input.cb)

sample_uid = snakemake.wildcards.sample_uid
min_genes_per_CB = int(snakemake.params.min_genes)
min_counts_per_CB = int(snakemake.params.min_counts)

# read 10X Genomics raw file
adata = sc.read_10x_h5(cellranger[1])
adata_cb = sc.read_10x_h5(cellbender)

adata.obs['sample_uid'] = sample_uid
adata.layers['background_removed'] = adata_cb.X

# filter out really low droplets
if snakemake.params.filter_cells:
    print("filtering cells based on {} genes and {} counts".format(min_genes_per_CB, min_counts_per_CB))
    sc.pp.filter_cells(adata, min_genes = min_genes_per_CB)
    sc.pp.filter_cells(adata, min_counts = min_counts_per_CB)

adata.write_h5ad(str(snakemake.output))
print("Done")




