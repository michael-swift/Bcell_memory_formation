import scanpy as sc
import anndata as ad
# get files from snake make
files = snakemake.input
# params from snakemake
min_genes = int(snakemake.params.min_genes)
min_counts = int(snakemake.params.min_counts)
filter_cells = snakemake.params.filter_cells
def load_and_filter(filename, min_genes=min_genes, min_counts=min_counts):
    adata = sc.read_h5ad(filename)
    adata.obs_names_make_unique()
    # filter very leniently
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
adata.obs_names_make_unique()
print(str(snakemake.output))
adata.write_h5ad(str(snakemake.output), compression = 'gzip')
print("Done!!")
