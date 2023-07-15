import scanpy as sc
import scanpy.external as sce
import scrublet as scr
import gc
import pathlib
import seaborn as sns

h5ad = str(snakemake.input)
output_file = str(snakemake.output.h5ad)
sample_id = snakemake.wildcards.sample_uid

def run_scrublet(adata, layers):
    for layer in layers:
        print(f"scrublet for {layer}")
        scrub = scr.Scrublet(adata.to_df(layer=layer))
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        # arbitrarily make take the top ventile as doublet associated (getting pretty wildly low doublet rates) reported
        adata.obs.loc[:, f"doublet_scores_{layer}"] = doublet_scores
        adata.obs.loc[:, f"predicted_doublets_{layer}"] = predicted_doublets
        adata.obs.loc[:, f"predicted_doublets_{layer}"] = adata.obs[f"predicted_doublets_{layer}"].astype(str)
    return adata

adata = sc.read_h5ad(h5ad)
adata.var_names_make_unique()
adata = run_scrublet(adata, layers = ['counts', 'cellbender_counts'])
adata.write_h5ad(output_file)
