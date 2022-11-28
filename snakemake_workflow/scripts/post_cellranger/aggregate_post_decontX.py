import scanpy as sc
import anndata as ad

files = glob.glob('*')

adatas = []
for i in files:
    adatas.append(sc.read_h5ad(i))

adata = ad.conat(adatas)

adata.write_h5ad("aggregated_decontX.h5ad")
