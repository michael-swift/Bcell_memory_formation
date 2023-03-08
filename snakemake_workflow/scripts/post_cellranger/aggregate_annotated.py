# reaggregate the tissue adata objects
obs_dfs = snakemake.input.obs_dfs
h5ads = snakemake.input.h5ads

for h5ad in h5ads:
    adata = sc.read_h5ad(h5ad)
    adatas.append(adata)
adata = ad.concat(adatas)
del adatas
gc.collect()
adata.write_h5ad(snakemake.output.full_h5ad, compression = 'gzip')

for obs_df in obs_dfs:
    df = pd.read_table(obs_df, index_col=0)
    dfs.append(df)
df = pd.concat(dfs)
df.to_csv(snakemake.output.gex_labels, sep = '\t')
