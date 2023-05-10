import scanpy as sc
import pandas as pd

cellranger = str(snakemake.input.cr)
cellbender = str(snakemake.input.cb)
sample_uid = snakemake.wildcards.sample_uid
min_genes_per_CB = int(snakemake.params.min_genes)
min_counts_per_CB = int(snakemake.params.min_counts)
samplesheets = snakemake.params.samplesheets
# read 10X Genomics raw file
print(cellranger)
adata = sc.read_10x_h5(cellranger)
adata_cb = sc.read_10x_h5(cellbender)

adata.obs["sample_uid"] = sample_uid
adata.layers["cellbender_counts"] = adata_cb.X
adata.layers["counts"] = adata.X

# filter out really low droplets
if snakemake.params.filter_cells:
    print(
        "filtering cells based on {} genes and {} counts".format(
            min_genes_per_CB, min_counts_per_CB
        )
    )
    sc.pp.filter_cells(adata, min_genes=min_genes_per_CB)
    sc.pp.filter_cells(adata, min_counts=min_counts_per_CB)


def add_samplesheet(samplesheets, adata):
    samplesheets.set_index("sample_uid", inplace=True)
    print(samplesheets.head)
    samplesheets = samplesheets[samplesheets.lib_type == "gex"]
    for column in samplesheets.columns:
        _dict = samplesheets[column].to_dict()
        adata.obs[column] = adata.obs.sample_uid.map(_dict)
        adata.obs.loc[:, column] = adata.obs[column].astype(str)
    return adata


samplesheets = pd.concat(
    [
        pd.read_table(
            "{}".format(x), dtype="str", sep="\t", index_col=0, engine="python"
        )
        for x in snakemake.params.samplesheets
    ],
    ignore_index=True,
)
print(samplesheets)
print(samplesheets["donor"])
print(samplesheets["sample_uid"])
adata = add_samplesheet(samplesheets, adata)
adata.obs = adata.obs[["donor", "tissue", "sample_uid"]]
adata.obs.to_csv("~/test.df.csv")
adata.write_h5ad(str(snakemake.output))
print("Done")
