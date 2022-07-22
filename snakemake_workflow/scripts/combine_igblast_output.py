import sys
import argparse
import pandas as pd
import numpy as np

###########################################################################

parser = argparse.ArgumentParser()
parser.add_argument('input_paths', nargs='+', 
                    help="list of files to combine")
parser.add_argument('-name', required="True", 
                    help="name for combined sample set")
parser.add_argument('-samplesheet', nargs='+', required=True, 
                    help="samplesheet containing sample metadata")
parser.add_argument('-outdir', default='.')

parser.add_argument('-pipeline_type', default='Mikey')
args = parser.parse_args()

filenames = args.input_paths
path_to_samplesheet = args.samplesheet
outdir = args.outdir
name = args.name
pipeline = args.pipeline_type
samplesheet = pd.concat([pd.read_table(
    "{}".format(x), dtype='str', sep='    ', engine='python')
    for x in path_to_samplesheet
], ignore_index=True)

# This is now hacked to make it work w my pipeline, where sample_id = sample_index = sample_uid
if pipeline == 'Mikey':
    samplesheet['sample_uid'] = samplesheet['sample_index']
else:
    samplesheet['sample_uid'] = samplesheet['library_uid'] + "-" + samplesheet['sample_index']
samplesheet.set_index('sample_uid', inplace=True)

metadata_columns = samplesheet.columns.values

all_data = []

for filename in filenames:
    sample_index = filename.split("/")[-1].split("_igblast.tsv")[0]
    sample_metadata = samplesheet.loc[sample_index]

    sample_df=pd.read_table(filename)
    for metadata_col in metadata_columns:
        sample_df[metadata_col] = sample_metadata[metadata_col]
    sample_df['sample_index'] = sample_index

    all_data.append(sample_df)

print("concatenating...")
df = pd.concat(all_data, ignore_index = True)
print('saving...')
df.to_csv("{path}/{file_id}_combined.tsv.gz".format(path=outdir, file_id=name), sep='\t', index=False)

