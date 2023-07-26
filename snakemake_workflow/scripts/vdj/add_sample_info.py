import argparse
import pandas as pd

###########################################################################

parser = argparse.ArgumentParser()
parser.add_argument('input', 
                    help="input file")
parser.add_argument('-samplesheet', nargs='+', required=True, 
                    help="samplesheet containing sample metadata")
parser.add_argument('-output', required=True)

args = parser.parse_args()

input_filename = args.input
path_to_samplesheet = args.samplesheet
output = args.output

samplesheet = pd.concat([pd.read_table(
    "{}".format(x), dtype='str', sep='\t', engine='python')
    for x in path_to_samplesheet
], ignore_index=True)

samplesheet = samplesheet.set_index('sample_uid')
samplesheet = samplesheet[samplesheet.lib_type=='vdj']

add_cols = ['donor', 
            'tissue', 
            'subanatomical_location',
            'sample_type',
            'ht_sample']

df = pd.read_table(input_filename)

# drop any samples removed from samplesheet
sample_in_samplesheet = df.sample_uid.isin(samplesheet.index)
if sample_in_samplesheet.sum() < df.shape[0]:
    print('Dropping cells from following samples:')
    print(df.sample_uid[~sample_in_samplesheet].unique())
df = df[sample_in_samplesheet]

for col in add_cols:
    suid_dict = samplesheet[col].to_dict()
    df[col] = df.sample_uid.map(suid_dict)
df.to_csv(output, sep='\t', index=False)


