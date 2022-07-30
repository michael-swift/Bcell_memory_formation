import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import argparse


parser = argparse.ArgumentParser(description='Combines VDJ and cell calls')
parser.add_argument('-cell_calls', 
                    help='cell call in tsv format')
parser.add_argument('-vdj_info', 
                    required=True)
parser.add_argument('-outdir', 
                    default='.', 
                    help="output directory (default: working directory)")
parser.add_argument('-outname', default='.',
                    help="name of output file. Default uses input file prefix up to first underscore.")

args = parser.parse_args()
CELL_CALLS = args.cell_calls
VDJ_INFO = args.vdj_info
OUTDIR = args.outdir
OUTNAME=args.outname

if OUTNAME == ".":
    OUTNAME = CELL_CALLS.split("/")[-1].split("_")[0]
OUTNAME += "_called_cells_" + "_".join(VDJ_INFO.split("/")[-1].split("_")[1:])

called_cell_df = pd.read_table(CELL_CALLS).rename(columns={'vdj_call': 'vdj_sequence'})
amplicon_df = pd.read_table(VDJ_INFO)

#drop columns that may be different for transcripts with the same vdj sequence
amplicon_df = amplicon_df.drop(['sequence_id',
                                'sequence_uid',
                                'rev_comp', 
                                'sequence', 
                                'SEQ_ID', 
                                'total_umis', 
                                'tissue', 
                                'sample_id', 
                                'sample_index',
                                'library_uid'], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.startswith("Unnamed")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.startswith("c_")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_start")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_end")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_cigar")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_support")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_btop")], axis=1)
amplicon_df = amplicon_df.drop_duplicates(ignore_index=True)
amplicon_df['j_call'] = amplicon_df['j_call'].map(lambda x : x.split('*')[0])
amplicon_df = amplicon_df.drop_duplicates(ignore_index=True)
called_cell_df = called_cell_df.merge(amplicon_df, on='vdj_sequence', how='left')

called_cell_df.to_csv('{}/{}'.format(OUTDIR, OUTNAME), index=False, sep='\t')
