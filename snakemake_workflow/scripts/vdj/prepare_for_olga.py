import sys
import argparse

import numpy as np
import pandas as pd

import itertools

######################## PATH CONFIG ################################

parser = argparse.ArgumentParser()
parser.add_argument('-input', required=True,
                    help="path to input data frame")
parser.add_argument('-outdir', required="True",
                    help="path to output directory")


######################## PARSE ARGS ################################

args = parser.parse_args()

input_path = args.input
outdir = args.outdir

######################################################################


full_df = pd.read_table(input_path, 
                        low_memory=False, 
                        usecols = ['vdj_sequence',
                                   'cdr3',
                                   'cdr3_aa',
                                   'fwr3',
                                   'fwr3_aa',
                                   'fwr4_aa',
                                   'fwr4',
                                   'donor', 
                                   'tissue',
                                   'c_call',
                                   'v_db_call',
                                   'j_call',
                                   'v_mismatch',
                                   'lineage_id'
                                   ])
full_df = full_df[full_df.vdj_sequence.notna()]

# find position fo conserved cysteine
full_df['C_pos'] = full_df.fwr3_aa.map(lambda x: len(x.split("C")[-1])+1)
# add conserved residues to cdr3
full_df['cdr3_full'] = full_df.apply(lambda x:x.fwr3[-3*(x.C_pos):],axis=1) + full_df['cdr3'] + full_df['fwr4'].map(lambda x:x[:3])
full_df['cdr3_aa_full'] = full_df.apply(lambda x:x.fwr3_aa[-(x.C_pos):],axis=1) + full_df['cdr3_aa'] + full_df['fwr4_aa'].map(lambda x:x[:1])


single_cdr3_lineages = full_df[(full_df.C_pos<=1) & (full_df.v_mismatch==0)].groupby(['donor','lineage_id'])['cdr3_aa_full'].nunique() ==1 
single_cdr3_lineages = single_cdr3_lineages.to_dict()
in_single_cdr3_lineage = full_df[(full_df.C_pos<=1) & (full_df.v_mismatch==0)].apply(lambda x: single_cdr3_lineages[(x.donor, x.lineage_id)], axis=1)
cdr3_df = full_df[(full_df.C_pos<=1) & (full_df.v_mismatch==0)][in_single_cdr3_lineage].groupby(['cdr3_full','cdr3_aa_full'])['donor'].unique()
cdr3_df = cdr3_df.reset_index()[['cdr3_full', 'cdr3_aa_full']]

list_df = np.array_split(cdr3_df, 200)
for i, x in enumerate(list_df):
    filename = "cdr3nt_chunk-" + str(i+1).zfill(3) + ".tsv"
    x.to_csv(f'{outdir}/{filename}', sep='\t', header=False)