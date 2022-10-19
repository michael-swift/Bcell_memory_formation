import scanpy as sc
from anndata import AnnData as ad
files = snakemake.input

for i, h5 in enumerate(files):
    if i==0:
        adata = sc.read_10x_h5(h5)
        # parse file name for sample_uid
        sample_uid=h5.split("/")[-3]
        print(sample_uid)
        adata.obs["sample_uid"] = sample_uid
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        # filter very leniently
        sc.pp.filter_cells(adata, min_genes = 50)
        sc.pp.filter_cells(adata, min_counts = 200)
    else:
        sub_adata = sc.read_10x_h5(h5)
        # parse file name for sample_uid
        sample_uid=h5.split("/")[-3]
        print(sample_uid)
        sub_adata.obs["sample_uid"] = sample_uid
        sub_adata.obs_names_make_unique()
        sub_adata.var_names_make_unique()
        # filter very leniently
        sc.pp.filter_cells(sub_adata, min_genes = 50)
        sc.pp.filter_cells(sub_adata, min_counts = 200)
        adata = adata.concatenate(sub_adata)
print(str(snakemake.output))
adata.write_h5ad(str(snakemake.output))

    
