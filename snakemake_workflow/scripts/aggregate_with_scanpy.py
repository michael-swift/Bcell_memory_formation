import scanpy as sc
from anndata import AnnData as ad
files = snakemake.input
min_genes = snakemake.params.min_genes
min_counts = snakemake.params.min_counts

def load_and_filter(filename, min_genes=min_genes, min_counts=min_counts):
    adata = sc.read_10x_h5(filename)
    # parse file name for sample_uid
    sample_uid=h5.split("/")[-3]
    print(sample_uid)
    adata.obs["sample_uid"] = sample_uid
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    # filter very leniently
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_cells(adata, min_counts = min_counts)
    return adata

# construct aggregated object
for i, h5 in enumerate(files):
    if i==0:
        adata = load_and_filter(h5)
    else:
        sub_adata = load_and_filter(h5)
        adata = adata.concatenate(sub_adata)

print(str(snakemake.output))
adata.write_h5ad(str(snakemake.output))
print("Done!!")
return
