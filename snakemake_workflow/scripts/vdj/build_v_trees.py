import sys, os

from io import StringIO, BytesIO
from Bio import SeqIO

from subprocess import Popen, PIPE

import argparse

import time

import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='Wrapper for phylogenetic reconstruction')
parser.add_argument('input_file', help="tabular file containing v MSAs")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-scratchdir', default = '.', help="scratch space for temp files (default: same as outdir)")
parser.add_argument('-samplename', default = 'default')
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

tmp_file = "{}/temp_file.fasta".format(scratchdir)

#########################################################################################################

v_seq_df = pd.read_table(FILENAME, usecols=['lineage_id', 'v_seq_id', 'v_sequence_msa'])
v_seq_df = v_seq_df.drop_duplicates(ignore_index=True)

#now pipe to fasttree
v_seq_df['fasta_records'] = ">" + v_seq_df.v_seq_id.astype(str) \
					 + "\n" + v_seq_df.v_sequence_msa

start = time.time()
sys.stderr.write("Piping sequences to FastTree...\n")

tree_list = []
last_time = 0
num_lineages = len(v_seq_df.lineage_id.unique())
for it, lineage in enumerate(v_seq_df.lineage_id.value_counts().index):
	subset = v_seq_df[v_seq_df.lineage_id == lineage]

	current_time = int(time.time() - start)
	if (len(subset) > 1000 ) or ((current_time % 10 == 0) and (current_time != last_time)):
		sys.stderr.write("[{}s] Inferring phylogeny for sequences in lineage {}/{} (n={})\n".format(int(time.time() - start), 
																				  it, num_lineages,
																				  len(subset)))
		last_time = current_time

	query_string = "\n".join(subset.fasta_records)

	with open("{}".format(tmp_file), 'w') as temp_buffer:
		temp_buffer.write(query_string)

	cmd = 'FastTree -nt -gtr {} '.format(tmp_file)
	p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
	out, err = p.communicate()

	out_nwk = str(out, encoding='utf-8')
	fasttree_log = str(err, encoding='utf-8')

	tree_list.append((lineage, out_nwk))
# # delete temporary file
os.remove(tmp_file)

with open("{}/{}_v_trees.tsv".format(outdir, samplename), 'w') as writer:
	writer.write("lineage_id\tv_phylogeny\n")
	for l, nwk in tree_list:
		writer.write("{}\t{}\n".format(l, nwk))
# v_seq_df.to_csv('{}/{}_vtrees.tsv.gz'.format(outdir,samplename), sep = '\t', index = False)
sys.stderr.write("Output written to: {}/{}_v_trees.tsv \n".format(outdir, samplename))
