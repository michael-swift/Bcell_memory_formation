import numpy as np
import pandas as pd
from copy import deepcopy

old = pd.read_table('old_samplesheet.csv', sep=',')
old = old.set_index(['sample_uid', 'lib_type'])
old['sample_index'] = 'SI-TT-' + old['sample_index']

rename_table = pd.read_table('rename_sample_uids.tsv')


new = deepcopy(rename_table)
new = new.rename(columns={'lib type':'lib_type'})
new['lib_type'] = new['lib_type'].str.lower()
new = new.set_index(['sample_uid', 'lib_type'])

new = new.merge(old, left_index=True, right_index=True, how='left').reset_index()

new['donor'] = new.sample_uid_final.map(lambda x: x.split("_")[0])
new['sample_type'] = new.sample_uid_final.map(lambda x: x.split("_")[1])
new['subanatomical_location'] = new.sample_uid_final.map(lambda x: x.split("_")[2].rstrip('over'))
tissue_dict = {'SDLN1':'LN', 
               'SDLN2':'LN', 
               'SDLN3':'LN',
               'MELN1':'LN',
               'MELN':'LN', 
               'BM-HL':'HL',
               'BM-HR':'HR'}
new['tissue'] = new.subanatomical_location.map(lambda x: tissue_dict.get(x,x))

new['ht_sample'] = new.sample_uid_final.str.contains('HT')
new['HT'] = new['HT'].astype(bool)


new_sample_info = pd.read_table('new_sample_info.tsv').set_index(['Sample', 'lib_type'])
new_sample_info['sample_index'] = 'SI-TT-' + new_sample_info['sample_index']

cell_count_dict = new_sample_info['expected_cells'].to_dict()
sample_index = new_sample_info['sample_index'].to_dict()

new.loc[new.expected_cells_thousands.isna(),
        'expected_cells_thousands'] = new.apply(lambda x: cell_count_dict.get((x.sample_uid,x.lib_type)),
                                                                                     axis=1)
new.loc[new.sample_index.isna(),
        'sample_index'] = new.apply(lambda x: sample_index.get((x.sample_uid,x.lib_type)),
                                                                                     axis=1)

new['expected_cells_thousands'] = new['expected_cells_thousands'].astype(int)
new = new[['sample_uid_final', 
            'lib_type', 
            'sample_index', 
            'expected_cells_thousands', 
            'donor', 
            'sample_type',
            'subanatomical_location', 
            'tissue', 
            'ht_sample']]
new = new.rename(columns={'sample_uid_final':'sample_uid'})
new = new.sort_values(['sample_uid', 'lib_type'])
new.to_csv('final_samplesheet.tsv', sep='\t', index=False)
