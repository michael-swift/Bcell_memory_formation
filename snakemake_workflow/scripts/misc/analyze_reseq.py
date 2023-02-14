import glob
import pandas as pd

# start from swiftShare dir
files = glob.glob("tabula_bursa_*pipeline*/per_sample/cellranger/*/outs/metrics_summary.csv")
dfs = []
for f in files:
    # parse file name
    split = f.split("/")
    if "_0" in split[0]:
        pipeline = "shallow"
    else:
        pipeline = "resequenced"
    sample = split[3]
    df = pd.read_csv(f, thousands = r',', sep = ',')
    df['pipeline'] = pipeline
    df['sample_uid'] = sample
    dfs.append(df)
df = pd.concat(dfs)
df[df.pipeline == 'resequenced'].sample_uid
resequenced = df[df.pipeline == 'resequenced'].sample_uid
df[df.sample_uid.isin(resequenced)]
sub_df = df[df.sample_uid.isin(reseqd)]
sub_df.groupby('pipeline')["Estimated Number of Cells"].mean()
