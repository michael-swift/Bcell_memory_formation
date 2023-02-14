import sys
from functools import partial

from Bio import SeqIO, bgzf
import regex

import gzip
import argparse

import time

import pandas as pd
import numpy as np

from pacbio_vdj_utils.cluster_vdj import *

import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description='Cluster lineages simple')
parser.add_argument('input_file', help="airr-formatted file")
parser.add_argument('-outdir', default='.', help="output directory (default: working directory)")
parser.add_argument('-samplename', default='.', help='samplename (default: filename retained)')
parser.add_argument('-max_distance_from_vstart', type=int, required=False, default=0)
parser.add_argument('-max_distance_from_cstart', type=int, required=False, default=0)
parser.add_argument('-fractional_cutoff', type=float, required=False, default=0.15)

args = parser.parse_args()

FILENAME = args.input_file
outdir = args.outdir
v_sstart_max = args.max_distance_from_vstart
c_sstart_max = args.max_distance_from_cstart
fractional_cutoff = args.fractional_cutoff
samplename = args.samplename

if samplename==".":
    samplename=FILENAME.split("/")[-1].split('.')[0]

#########################################################################################################

nrows = None
# read only relevant columns first

df = pd.read_table(FILENAME, usecols=['sequence', 
                              'v_family',
                              'v_sequence_alignment',
                              'j_sequence_alignment',
                              'v_sequence_start', 
                              'v_germline_start',
                              'v_germline_alignment',
                              'j_germline_alignment',
                              'cdr3_start', 
                              'cdr3_end', 
                              'cdr3',
                              'j_sequence_end'], nrows=nrows)

df['cdr3_length'] = df['cdr3'].str.len()

sys.stderr.write("Starting with %d reads...\n" % df.shape[0])

df = df[df.v_germline_start <= v_sstart_max + 1]
try:
    df = df[df['c_sstart'] <= c_sstart_max + 1]
except KeyError:
    pass
sys.stderr.write("Keeping {} reads that contain entire VDJ region...\n".format(df.shape[0]))

#df = df[df['v_family'].str.startswith("IGH")]
sys.stderr.write("Subsetting to {} reads that map to heavy chain...\n".format(df.shape[0]))

df['vdj_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start)-1:int(x.j_sequence_end)], axis = 1)

df['v_templated_len'] = df['cdr3_start'] - df['v_sequence_start']
df['j_templated_len'] = df['j_sequence_end'] - df['cdr3_end']

# drop reads with gap in templated alignment
df['sum_gap_len'] = df.v_sequence_alignment.map(lambda x: x.count("-")) + df.j_sequence_alignment.map(lambda x: x.count("-"))
df = df[df['sum_gap_len'] < 1]

sys.stderr.write("Keeping {} reads with no deletions in templated sequences...\n".format(df.shape[0]))


df['sum_insertion_len'] = df.v_germline_alignment.map(lambda x: x.count("-")) + df.j_germline_alignment.map(lambda x: x.count("-"))
df = df[df['sum_insertion_len'] < 1]

sys.stderr.write("Keeping {} reads with no insertions in templated sequences...\n".format(df.shape[0]))

unique_vdjs = df[['v_family',
                  'vdj_sequence', 
                  'v_templated_len', 
                  'j_templated_len',
                  'cdr3', 
                  'cdr3_length']].drop_duplicates(ignore_index=True)
sys.stderr.write("Clustering {} unique variable sequences...\n".format(unique_vdjs.shape[0]))

total_clusters = 0
unique_vdjs['cluster_id'] = -1

start = time.time()
distance_to_nearest = []
for vfam in unique_vdjs.v_family.value_counts().index:
	for cdr3_len in unique_vdjs[unique_vdjs['v_family'] == vfam].cdr3_length.value_counts().index:
		
		subset = unique_vdjs[(unique_vdjs['v_family'] == vfam) & (unique_vdjs['cdr3_length'] == cdr3_len)]
		sys.stderr.write("Current subset: v_family={}, cdr3_length={}, n={}\n".format(vfam, cdr3_len, len(subset)))

		cdr3_seqs = subset.cdr3.values
		subset_idx = subset.index

		Ds = get_hamming_distances(cdr3_seqs, squareform=True)/cdr3_len

		n = Ds.shape[0]
		diagonal_indices = np.diag_indices(n)
		Ds[diagonal_indices] = 1
		closest_d = np.min(Ds, axis=1)
		distance_to_nearest.extend(list(closest_d))
		Ds[diagonal_indices] = 0

		upper_triangular_indices = np.triu_indices(n, 1)
		condensed_distance_matrix = Ds[upper_triangular_indices]

		new_cluster_ids = get_cluster_ids(condensed_distance_matrix, 
										 cutoff=fractional_cutoff, 
										 method='single')
		unique_vdjs.loc[subset_idx, 'cluster_id'] = total_clusters + new_cluster_ids
		total_clusters += max(new_cluster_ids) + 1

		sys.stderr.write("[%ds]  found %d new cdr3 clusters, continuing...\n" % (time.time() - start, new_cluster_ids.max() + 1))

sys.stderr.write("clustering templated sequences within cdr3 clusters...\n")

unique_vdjs['templated_v'] = unique_vdjs.apply(lambda x: x.vdj_sequence[0:int(x.v_templated_len)], axis=1)
unique_vdjs['templated_j'] = unique_vdjs.apply(lambda x: x.vdj_sequence[-int(x.j_templated_len):], axis=1)
total_lineages = 0
unique_vdjs['lineage_ids'] = -1

for cluster in unique_vdjs.cluster_id.value_counts().index:
    subset = unique_vdjs[unique_vdjs.cluster_id == cluster]
    sys.stderr.write("Current subset: cluster_id={}, n={}\n".format(cluster, len(subset)))

    templ_v_seqs = subset.templated_v.values
    templ_j_seqs = subset.templated_j.values
    subset_idx = subset.index
    n = len(subset.index)

    Ds_V = get_pairwise_distances(templ_v_seqs, squareform=True)
    Ds_J = get_pairwise_distances(templ_j_seqs, squareform=True)

    Ds = Ds_V + Ds_J

    lengths = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1):

            lengths[i,j] = max(len(templ_v_seqs[i] + templ_j_seqs[i]), len(templ_v_seqs[j] + templ_j_seqs[j]))
            lengths[j,i] = lengths[i,j]

    Ds = Ds/lengths

    upper_triangular_indices = np.triu_indices(n, 1)
    condensed_distance_matrix = Ds[upper_triangular_indices]
    new_lineage_ids = get_cluster_ids(condensed_distance_matrix, 
                                     cutoff=fractional_cutoff, 
                                     method='single')

    unique_vdjs.loc[subset_idx, 'lineage_ids'] = total_lineages + new_lineage_ids
    total_lineages += max(new_lineage_ids) + 1

    sys.stderr.write("[%ds]  found %d total lineages, continuing...\n" % (time.time() - start, new_lineage_ids.max() + 1))


fig, ax = plt.subplots()
bins = np.linspace(0,1,40)
counts, bins = np.histogram(distance_to_nearest, bins = bins)
bin_centers = (bins[1:] + bins[:-1])/2.
ax.plot(bin_centers, counts)
ax.set_yscale('log')
ax.axvline(x=fractional_cutoff, color='r')
ax.set_xlabel('Distance to nearest cdr3 (fractional)')
ax.set_ylabel('Number of sequences')
fig.savefig('{}/{}_distance_to_nearest_cdr3.pdf'.format(outdir, samplename), bbox_inches = 'tight')
plt.close(fig)


vdj_families = unique_vdjs[['vdj_sequence','lineage_ids']].set_index('vdj_sequence').to_dict()
vdj_families = vdj_families['lineage_ids']

# read whole dataframe
df = pd.read_table(FILENAME, nrows = nrows)
df['cdr3_length'] = df['cdr3'].str.len()

df = df[df.v_germline_start <= v_sstart_max + 1]
try:
    df = df[df['c_sstart'] <= c_sstart_max + 1]
except KeyError:
    pass
df = df[df['v_family'].str.startswith("IGH")]
df['sum_gap_len'] = df.v_sequence_alignment.map(lambda x: x.count("-")) + df.j_sequence_alignment.map(lambda x: x.count("-"))
df = df[df['sum_gap_len'] < 1]
df['sum_insertion_len'] = df.v_germline_alignment.map(lambda x: x.count("-")) + df.j_germline_alignment.map(lambda x: x.count("-"))
df = df[df['sum_insertion_len'] < 1]

df['vdj_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start)-1:int(x.j_sequence_end)], axis = 1)
df['lineage_id'] = df.vdj_sequence.map(lambda x: vdj_families[x])


df.to_csv('{}/{}.tsv.gz'.format(outdir,samplename), sep = '\t', index = False)
