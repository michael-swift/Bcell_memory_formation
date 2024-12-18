import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Combines VDJ and cell calls')
parser.add_argument('-cell_calls', 
                    help='cell call in tsv format')
parser.add_argument('-vdj_info', 
                    required=True)
parser.add_argument('-outdir', 
                    default='.', 
                    help="output directory (default: working directory)")
parser.add_argument('-outname', default='.',
                    help="name of output file. Default uses input file prefix.")

args = parser.parse_args()
CELL_CALLS = args.cell_calls
VDJ_INFO = args.vdj_info
OUTDIR = args.outdir
OUTNAME=args.outname

if OUTNAME == ".":
    OUTNAME = CELL_CALLS.split("/")[-1].split(".")[0]


called_cell_df = pd.read_table(CELL_CALLS).rename(columns={'vdj_call': 'vdj_sequence'})
amplicon_df = pd.read_table(VDJ_INFO)

#drop columns that may be different for transcripts with the same vdj sequence
amplicon_df = amplicon_df.drop(['sequence_id',
                                'sequence_uid',
                                'reads',
                                'rev_comp', 
                                'sequence', 
                                'umis',
                                'lib_type', 
                                'sample_index',
                                'expected_cells_thousands',
                                'donor',
                                'sample_type',
                                'subanatomical_location',
                                'tissue',
                                'ht_sample'], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.startswith("Unnamed")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.startswith("c_")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_start")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_end")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_cigar")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_support")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_btop")], axis=1)

# drop 10X barcode and contig_id columns because their nomenclature is inconsistent with our cell call dataframe
amplicon_df = amplicon_df.drop(['barcode', 'contig_id'], axis=1)
# drop 10X annotations, since these may be different for the same amplicon
tenX_annotations=[x for x in amplicon_df.columns if (x.endswith("10X") or x.startswith("10X"))]
amplicon_df = amplicon_df.drop(tenX_annotations, axis=1)
amplicon_df = amplicon_df.drop_duplicates(ignore_index=True)

amplicon_df['j_call'] = amplicon_df['j_call'].map(lambda x : x.split('*')[0])
amplicon_df = amplicon_df.drop_duplicates(ignore_index=True)

# merge and write
called_cell_df = called_cell_df.merge(amplicon_df, on=['sample_uid', 'vdj_sequence', 'locus'], how='left')

called_cell_df.to_csv(f'{OUTDIR}/{OUTNAME}.tsv.gz', index=False, sep='\t')
