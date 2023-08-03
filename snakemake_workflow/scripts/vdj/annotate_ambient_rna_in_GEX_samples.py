import sys
import argparse

import numpy as np
import pandas as pd

import itertools

######################## PATH CONFIG ################################

parser = argparse.ArgumentParser()
parser.add_argument('-input', required=True,
                    help="path to input data frame")
parser.add_argument('-output', required="True", 
                    help="path to output data frame")
parser.add_argument('-min_umis', default=500, type=int)
parser.add_argument('-min_frac', default=0.25, type=float)


######################## PARSE ARGS ################################

args = parser.parse_args()

input_path = args.input
output_path = args.output
MIN_UMIS = args.min_umis
FRAC_TOTAL = args.min_frac
######################## ############ ################################


df = pd.read_table(input_path, low_memory=False)

################################################################################################################

df['has_vdj'] = df.vdj_sequence.notna()

is_ambient_vdj = df[(df.has_vdj)
                     & (df.sample_uid_gex.notna())].groupby(['sample_uid_gex', 
                                                             'vdj_sequence'])['probable_hq_single_not_b_cell'].sum()
is_ambient_vdj = is_ambient_vdj > 0
is_ambient_vdj = is_ambient_vdj.to_dict()

df['vdj_in_ambient'] = False
df.loc[(df.has_vdj)
       & (df.sample_uid_gex.notna()), 
       'vdj_in_ambient'] = df[(df.has_vdj)
                              & (df.sample_uid_gex.notna())].apply(lambda x: 
                                            is_ambient_vdj[(x.sample_uid_gex, 
                                                            x.vdj_sequence)], axis=1)
vdjs = df[df['vdj_in_ambient']].groupby(['sample_uid_gex', 
                                      'vdj_sequence'])['n_umis'].agg(list).reset_index()

vdjs['n_umis'] = vdjs['n_umis'].map(lambda x: np.asarray(x))
vdjs['total'] = vdjs['n_umis'].map(lambda x: sum(x))
vdjs['max'] = vdjs['n_umis'].map(lambda x: max(x))
vdjs = vdjs[vdjs['total'] > vdjs['max']]

ambient_vdj_dict = vdjs.set_index(['sample_uid_gex','vdj_sequence'])[['total', 'max']].to_dict()

drop_column_list = ['vdj_umi_total', 'vdj_umi_max']


df['vdj_umi_max'] = df.apply(
                             lambda x: ambient_vdj_dict['max'].get(
                                    (x.sample_uid_gex, x.vdj_sequence), -1),
                              axis=1)
df['vdj_umi_total'] = df.apply(
                               lambda x: ambient_vdj_dict['total'].get(
                                        (x.sample_uid_gex, x.vdj_sequence), -1),
                               axis=1)

source_filter = ((df.n_umis == df.vdj_umi_max).astype(bool) 
                  & (df.n_umis > MIN_UMIS).astype(bool)
                  & (df.n_umis > (FRAC_TOTAL * df.vdj_umi_total)).astype(bool))


df['is_ambient_source'] = False
df.loc[source_filter, 'is_ambient_source'] = True


source_dict = df[source_filter].groupby('sample_uid_gex')['vdj_sequence'].agg(list)
source_cell_dict = df[source_filter][['sample_uid_gex','vdj_sequence','cb']]
source_cell_dict = source_cell_dict.set_index(['sample_uid_gex', 'vdj_sequence'])
source_cell_dict = source_cell_dict['cb'].to_dict()

df['vdj_is_from_ambient'] = df.apply(lambda x: 
                        x.vdj_sequence in source_dict.get(x.sample_uid, []),
                                     axis=1
                                    )
df.loc[source_filter, 'vdj_is_from_ambient'] = False

df['source_cell'] = np.nan
df['source_cell'] = df[df.vdj_is_from_ambient].apply(lambda x: 
                    source_cell_dict[(x.sample_uid_gex, x.vdj_sequence)],
                                                    axis=1)
df.loc[df.sample_uid_gex.isna(), 'is_ambient_source'] = np.nan
df.loc[df.sample_uid_gex.isna(), 'vdj_in_ambient']= np.nan

df = df.drop(drop_column_list, axis=1)
df.to_csv(output_path, sep='\t', index=False)

