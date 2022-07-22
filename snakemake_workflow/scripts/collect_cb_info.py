import numpy as np
import pandas as pd

import argparse

import time

#####################################################    

def add_cb_umi(cb_dict,cb_umi,count):
    cb, umi = cb_umi.split("_")
    cb_dict[cb] = cb_dict.get(cb,{})
    cb_dict[cb][umi] = cb_dict[cb].get(umi,0) + count

#####################################################    

parser = argparse.ArgumentParser(description='Parses tsv containing links between cell barcodes, umis, and seqids')
parser.add_argument('input_file')
parser.add_argument('-outdir', default='.', help="output directory (defualt: working directory)")
args = parser.parse_args()

FILENAME = args.input_file
outdir = args.outdir

sample_uid=FILENAME.split("/")[-1].split('_consensus-pass_seq_ids_barc-des_corrected.tsv')[0]

#####################################################

cb_df = pd.read_table(FILENAME)

cb_df['cell_barcodes'] = cb_df.SEQ_ID_INFO.map(lambda x: [y for y in x.split("|")[0].split("=")[-1].split(",")])
cb_df['reads'] = cb_df.SEQ_ID_INFO.map(lambda x: [int(y) for y in x.split("|")[-1].split("=")[-1].split(",")])

cb_seq_id_dict = {}

cb_umi_dict = {}

for it, row in cb_df.iterrows():
    for it2 in range(len(row.cell_barcodes)):
        add_cb_umi(cb_umi_dict,row.cell_barcodes[it2],row.reads[it2])
        cb = row.cell_barcodes[it2].split("_")[0]
        if cb in cb_seq_id_dict.keys():
            cb_seq_id_dict[cb] = cb_seq_id_dict[cb] + [row.SEQ_ID]
        else:
            cb_seq_id_dict[cb] = [row.SEQ_ID]

for k, v in cb_seq_id_dict.items():
    cb_seq_id_dict[k] = ",".join([str(x) for x in set(v)])

for k, v in cb_umi_dict.items():
    umi_data, read_data = zip(*[(uk,str(uv)) for uk, uv in v.items()])
    cb_umi_dict[k] = [",".join(umi_data), ",".join(read_data)]


cb_df = pd.DataFrame.from_dict(cb_seq_id_dict, orient='index').reset_index()
cb_df = cb_df.rename(columns={'index':'CB', 0:'SEQ_IDs'})

cb_df.to_csv('{}/{}_cb_seqids_whitelisted.tsv.gz'.format(outdir, sample_uid), sep='\t', index=False)


cb_df = pd.DataFrame.from_dict(cb_umi_dict, orient='index').reset_index()
cb_df = cb_df.rename(columns={'index':'CB', 0:'UMIs', 1:'READS'})

cb_df.to_csv('{}/{}_cb_umis_whitelisted.tsv.gz'.format(outdir, sample_uid), sep='\t', index=False)

    