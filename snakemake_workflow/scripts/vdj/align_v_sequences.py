import sys, os

from io import StringIO, BytesIO
from Bio import SeqIO

from subprocess import Popen, PIPE

import argparse

import time

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Wrapper for v sequence MSA')
parser.add_argument('input_file', help="airr-formatted file")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-scratchdir', default = '.', help="scratch space for temp files (default: same as outdir)")
parser.add_argument('-samplename', default = 'default')
parser.add_argument('-germline_db', default = 'path to germline database')
parser.add_argument('-threads', type=int, default=1)

args = parser.parse_args()

FILENAME = args.input_file
outdir = args.outdir
if args.scratchdir == '.':
	scratchdir = outdir
else:
	scratchdir = args.scratchdir

if args.samplename == 'default':
	samplename=FILENAME.split("/")[-1].split('.tsv.gz')[0]
else:
	samplename = args.samplename
germline_db_path = args.germline_db

THREADS = args.threads

tmp_in_file = "{}/in.fasta".format(scratchdir)
tmp_align_file = "{}/aln.fasta".format(scratchdir)
#########################################################################################################

target_seq = "AGGTGCAGTTGGTGGAATCTGGGGGAGGCGTGGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTCGCTATGCTATACACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAACAATCAATATTACGCAGACACCATGCAGGGCCGATTCACCATCTCCAGAGACAATTCCGAGAACATGCTGTATCTGCAAATGAACAGCCTGGGACCTGAGGACACAGGTCTTTATTATTGTG"
germline_sequences = SeqIO.to_dict(SeqIO.parse(germline_db_path, "fasta"))
germline_sequences = {k:str(v.seq) for k, v in germline_sequences.items()}

v_seq_df = pd.read_table(FILENAME, usecols=['locus','v_sequence', 'lineage_id', 'v_db_call'])
v_seq_df = v_seq_df.drop_duplicates(ignore_index=True)
v_seq_df = v_seq_df[v_seq_df.locus=='IGH']
#v_seq_df = v_seq_df.iloc[1000:1200,:]

all_v_seqs = v_seq_df.v_sequence.unique()
all_lineages = v_seq_df.lineage_id.unique()
print('started with {} unique data v seqs'.format(len(all_v_seqs)), file=sys.stderr)
print('these seqs belong to {} lineages'.format(len(all_lineages)), file=sys.stderr)

# add germline sequence corresponding to every gene to each lineage
lineage_majority_germline = v_seq_df.groupby('lineage_id')['v_db_call'].agg(lambda x: x.value_counts().index[0]).reset_index()
lineage_majority_germline['v_sequence'] = lineage_majority_germline.v_db_call.map(germline_sequences)
v_seq_df = pd.concat([v_seq_df, lineage_majority_germline])

# delete germline from lineages in which it existed already to avoid duplicates
v_seq_df = v_seq_df.drop_duplicates(ignore_index=True)

v_seq_ids = {seq: "V"+str(it) for it, seq in enumerate(v_seq_df.v_sequence.unique())}
print(v_seq_df.v_db_call.nunique(), "unique v db calls", file=sys.stderr)
print(v_seq_df.v_sequence.nunique(), "unique v sequences", file=sys.stderr)
print(v_seq_df.lineage_id.nunique(), "unique lineages", file=sys.stderr)
v_seq_df['v_seq_id'] = v_seq_df['v_sequence'].map(lambda x: v_seq_ids[x])
v_seq_df['v_seq_id'] = v_seq_df['lineage_id'].astype(str) + "_" + v_seq_df['v_seq_id']
#now pipe to muscle
v_seq_df['fasta_records'] = ">" + v_seq_df.v_seq_id.astype(str) \
					 + "\n" + v_seq_df.v_sequence

start = time.time()
sys.stderr.write("Piping sequences to muscle...\n")

aligned_seq_records = {}
last_time = 0
num_lineages = len(v_seq_df.lineage_id.unique())
for it, lineage in enumerate(v_seq_df.lineage_id.value_counts().index):
    subset = v_seq_df[v_seq_df.lineage_id == lineage]
    unique_seqs = subset.fasta_records.unique()

    current_time = int(time.time() - start)
    if (len(unique_seqs) > 1000 ) or ((current_time % 10 == 0) and (current_time != last_time)):
        sys.stderr.write("[{}s] Aligning sequences in lineage {}/{} (n={})\n".format(int(time.time() - start), 
	    																		  it, num_lineages,
																				  len(unique_seqs)))
        last_time = current_time

    query_string = "\n".join(unique_seqs)
     
    with open("{}".format(tmp_in_file), 'w') as temp_buffer:
        temp_buffer.write(query_string)
    cmd = 'muscle -align {} -output {} -nt -threads {}'.format(tmp_in_file, tmp_align_file, THREADS)
     
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
    out, err = p.communicate()
    out_fasta = str(out, encoding='utf-8')
    muscle_log = str(err, encoding='utf-8')

    aligned_seq_records.update({rec.id:str(rec.seq) for rec in SeqIO.parse(open(tmp_align_file,'r'), "fasta")})
# delete temporary files
#os.remove(tmp_in_file)
#os.remove(tmp_align_file)

v_seq_df['v_sequence_msa'] = v_seq_df['v_seq_id'].map(lambda x : aligned_seq_records.get(x,''))
v_seq_df = v_seq_df[['v_sequence', 'lineage_id', 'v_seq_id', 'v_sequence_msa']]

# almost ready to write, just annotate new germline sequences added to lineages
# in which raw data did not contain new germline
df = pd.read_table(FILENAME, usecols=['v_sequence', 'lineage_id']).drop_duplicates(ignore_index=True)
df['source'] = 'data'

v_seq_df = v_seq_df.merge(df, on=['v_sequence','lineage_id'], how='left')
v_seq_df['source'] = v_seq_df.source.fillna('inferred')
v_seq_df['v_seq_id'] = v_seq_df['v_seq_id'] + "_" + v_seq_df['source']
v_seq_df = v_seq_df[['v_sequence', 'lineage_id', 'v_seq_id', 'v_sequence_msa']]

found = np.asarray([seq in v_seq_df.v_sequence.values for seq in all_v_seqs], dtype=bool)
print("found {} of {} v sequences in final df".format(sum(found), len(all_v_seqs)), file=sys.stderr)
found_lin = np.asarray([lin in v_seq_df.lineage_id.values for lin in all_lineages], dtype=bool)
print("found {} of {} lineages in final df".format(sum(found_lin), len(all_lineages)), file=sys.stderr)
if found.sum() < len(found):
    print("first missing sequence ",all_v_seqs[~found][0], file=sys.stderr)
if found_lin.sum() < len(found_lin):
    print("first missing lineage ", all_lineages[~found_lin][0], file=sys.stderr)

v_seq_df.to_csv('{}/{}_vmsa.tsv.gz'.format(outdir,samplename), sep = '\t', index = False)
sys.stderr.write("Output written to: {}/{}_vmsa.tsv.gz \n".format(outdir, samplename))
