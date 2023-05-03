import scanpy as sc
import pandas as pd
import anndata as ad
import gc

# reaggregate the tissue adata objects
h5ads = snakemake.input.h5ads
print(len(h5ads), "duplication by snakemake causing a lot of h5ads")
h5ads = set(h5ads)
print(len(h5ads), "set of h5ads")
adatas = []
for h5ad in h5ads:
    print(h5ad)
    adata = sc.read_h5ad(h5ad)
    adatas.append(adata)
adata = ad.concat(adatas)
del adatas
gc.collect()
adata.obs_names_make_unique()
adata.write_h5ad(snakemake.output.full_h5ad, compression="gzip")
adata.obs.to_csv(snakemake.output.adata_obs)
