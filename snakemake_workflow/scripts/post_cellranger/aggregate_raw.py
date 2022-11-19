import scanpy as sc
import anndata as ad
# get files from snake make
files = snakemake.input
# params from snakemake
min_genes = int(snakemake.params.min_genes)
min_counts = int(snakemake.params.min_counts)

def load_and_filter(filename, min_genes=min_genes, min_counts=min_counts):
    adata = sc.read_h5ad(filename)
    # filter very leniently
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_cells(adata, min_counts = min_counts)
    return adata

adata_list = []

# construct aggregated object
for i, file_name in enumerate(files):
    adata_list.append(load_and_filter(file_name))
    print(i, file_name, " anndatas added to list")

adata_agg = ad.concat(adata_list)
#  filter genes which are largely undetected across all cells
sc.pp.filter_genes(adata_agg, min_counts = 10)
print(str(snakemake.output))
adata.write_h5ad(str(snakemake.output))
print("Done!!")
