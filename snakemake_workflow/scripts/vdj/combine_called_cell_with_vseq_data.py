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
                                'tissue',
                                'sample_index',
                                'donor',
                                'sample_type',
                                'sample_descriptor',
                                'tissue',
                                'BA CONC (ng/uL)',
                                'LIB QUALITY (0-3)',
                                'LIB CONC (ng/uL)',
                                'LIB MOLARITY (nM)',
                                'LIB AVG LEN (bp)',
                                'Pool weight',
                                'expect cells (k)',
                                'lib_type'], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.startswith("Unnamed")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.startswith("c_")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_start")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_end")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_cigar")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_support")], axis=1)
amplicon_df = amplicon_df.drop([x for x in amplicon_df.columns if x.endswith("_btop")], axis=1)

# drop 10X barcode and contig_id columns because their nomenclature is inconsistent with our cell call dataframe
amplicon_df = amplicon_df.drop(['barcode', 'contig_id'], axis=1)

amplicon_df = amplicon_df.drop_duplicates(ignore_index=True)
amplicon_df['j_call'] = amplicon_df['j_call'].map(lambda x : x.split('*')[0])
amplicon_df = amplicon_df.drop_duplicates(ignore_index=True)

# aggregate 10X annotations, since these may be different for the same amplicon
tenX_annotations=[x for x in amplicon_df.columns if x.endswith("10X")]
for x in tenX_annotations:
    df[x] = df[x].astype(str)

other_annotations = [x for x in amplicon_df.columns if not (x.endswith("10X"))]
amplicon_df = amplicon_df.groupby(other_annotations)[tenX_annotations].agg(",".join).reset_index()

# merge and write
called_cell_df = called_cell_df.merge(amplicon_df, on=['sample_uid', 'vdj_sequence', 'locus'], how='left')

called_cell_df.to_csv(f'{OUTDIR}/{OUTNAME}.tsv.gz', index=False, sep='\t')
