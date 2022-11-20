import scanpy as sc

adata = sc.read_h5ad(snakemake.input[0])
adata.obs['sample_uid'] = snakemake.wildcards.sample_uid
print("setting var index to feature names")
adata.var.set_index(adata.var['feature_name'], inplace = True)
adata.write_h5ad(snakemake.output[0])
print("done!")
