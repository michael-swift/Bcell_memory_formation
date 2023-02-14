import argparse

import sys
import time
from copy import deepcopy

import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from Levenshtein import distance

from pacbio_vdj_utils.trie_dict_utils import read_trie, search

#####################################################

parser = argparse.ArgumentParser(description='Calls cells based on VDJC sequence data emulsion based sc experiments')
parser.add_argument('input',
                    help='VDJC sequence file in tsv format')
parser.add_argument('-whitelist',
                    help="cell barcode whitelist")
parser.add_argument('-outdir',
                    default='.',
                    help="output directory (default: working directory)")
parser.add_argument('-max_VDJ_error',
                    default=10,
                    help="max error in vdj part of sequence")
parser.add_argument('-min_droplet_fraction', default=0.1, type=float,
                    help="minimal fraction of UMIs within droplet to call a cell")
parser.add_argument('-outname', default='.',
                    help="name of output file. Default uses input file prefix.")

args = parser.parse_args()
VDJC_SEQUENCE_TSV = args.input
WHITELIST_FILE=args.whitelist
MAX_ERROR_IN_VDJ = args.max_VDJ_error
OUTDIR = args.outdir
OUTNAME=args.outname
MIN_DROPLET_FRACTION = args.min_droplet_fraction

if OUTNAME == ".":
    OUTNAME=VDJC_SEQUENCE_TSV.split("/")[-1].split(".")[0]

#####################################################

WHITELIST = pd.read_table(WHITELIST_FILE, header=None)[0].values

start = time.time()
trie = read_trie(WHITELIST)
end = time.time()
elapsed = round((end - start)*1000)/1000
sys.stderr.write(f"Constructing whitelist trie took {elapsed}s\n")

#####################################################

df = pd.read_table(VDJC_SEQUENCE_TSV,
                    usecols=['barcode',
                             'locus',
                             'sequence_id',
                             'sample_uid',
                             'vdj_sequence',
                             'c_gene_10X',
                             'sequence',
                              'umis',
                              'reads'])

df['cb'] = df['barcode'].map(lambda x : x.split("-")[0])

sys.stderr.write('Retaining only whitelisted barcodes')

all_cell_barcodes = df['cb'].unique()

CBs = {k:None for k in all_cell_barcodes}

start = time.time()

EXACT = 0
CORRECTED = 0
TOTAL = 0

for cb in CBs.keys():
    TOTAL += 1
    if len(search(cb, trie, 0)) == 1:
        CBs[cb] = cb
        EXACT += 1
    else:
        single_edit = search(cb, trie, 1)
        if len(single_edit) == 1:
            CBs[cb] = single_edit[0][0]
            CORRECTED += 1

end = time.time()

elapsed = round((end - start)*1000)/1000
sys.stderr.write(f"Whitelisting and error correction took {elapsed}s\n")
sys.stderr.write(f"{TOTAL} barcodes processed\n")
sys.stderr.write(f"{(EXACT + CORRECTED)} barcodes whitelisted ({EXACT} exact)\n")
sys.stderr.write(f"{(TOTAL - EXACT - CORRECTED)} cell barcodes discarded\n")

df['cb_corrected'] = df['cb'].map(CBs)
df = df[df.cb_corrected.notna()]
df['cb'] = df['cb_corrected']

####################################################

sys.stderr.write('aggregating contigs belonging to same cell barcode...\n')

vdjc_df = df.groupby(['sample_uid',
                      'cb',
                      'locus',
                      'vdj_sequence',
                      'c_gene_10X'])[['umis', 'reads']].agg(sum)
cb_df = vdjc_df.reset_index().groupby(['sample_uid',
                                       'cb',
                                       'locus'])[['vdj_sequence',
                                                  'c_gene_10X',
                                                  'umis',
                                                  'reads']].agg(list)
cb_df = cb_df.reset_index()

for col in ['umis', 'reads']:
    cb_df[f'total_{col}'] = cb_df[col].map(lambda x: sum(x))
    cb_df[f'{col}'] = cb_df[col].apply(lambda x: np.array(x, dtype=int))
for col in ['vdj_sequence', 'c_gene_10X']:
    cb_df[col] = cb_df[col].apply(lambda x: np.array(x))

print('INFO: Dropping contigs with fewer than 2 UMIs', file=sys.stderr)

# cell barcode has to be associated with at least 2 UMIs
cb_df = cb_df[cb_df.total_umis > 1]

print('INFO: total number of valid contigs with at least 2 UMIs:', file=sys.stderr)
print(cb_df.groupby('sample_uid').size(), file=sys.stderr)

cb_df = cb_df.sort_values(by='total_umis', ascending=False)

########################################################
is_cell_dict = {}
is_multiplet_dict = {}
vdj_call_dict = {}
corrected_vdj_dict = {}

components = []

for it, row in cb_df.iterrows():

    unique_vdj_seqs = np.unique(np.asarray(row.vdj_sequence))
    umi_support = [(sum(row.umis[row.vdj_sequence == x]),x) for x in unique_vdj_seqs]
    umi_support = sorted(umi_support)

    umis_per_seq = np.asarray([x for x, y in umi_support])
    unique_vdj_seqs = [y for x, y in umi_support]

    max_umi_support = umi_support[-1][0]
    n = len(unique_vdj_seqs)
    D = np.zeros((n,n))

    corrected_vdj_dict.update({it:{}})
    # if multiple VDJ sequences, get distances
    if n > 1:
        for i in range(n):
            for j in range(i):
                d = distance(unique_vdj_seqs[i],unique_vdj_seqs[j])
                D[i,j] = d
                D[j,i] = d

        # first cut this into connected components
        # max distance between connected components is MAX_ERROR_IN_VDJ
        graph = csr_matrix(D < MAX_ERROR_IN_VDJ)
        n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)

        labels = np.asarray(labels)
        umi_support_per_component = np.zeros(n_components)
        for i in range(n_components):
            umi_support_per_component[i] = (umis_per_seq[labels == i]).sum()

        # keep only components supported by multiple UMIs, and no fewer than 10% of droplet content
        # if more than one connected component require at least 4 UMIs
        UMI_support_threshold = 1 + 2 * int(n_components > 1)
        fractional_umi_support = umi_support_per_component/umi_support_per_component.sum()
        good_components = np.arange(n_components)[(fractional_umi_support > MIN_DROPLET_FRACTION) \
                          & (umi_support_per_component > UMI_support_threshold)]
        components.append((len(good_components), n_components))

        if len(good_components) > 0:
            is_cell_dict[it] = True

            if len(good_components) > 1:
                is_multiplet_dict[it] = True
            else:
                is_multiplet_dict[it] = False

            cell_vdj_calls = []

            for good_component in good_components:
                good_component_locs = np.arange(len(labels))[labels == good_component]
                D_subset = D[np.ix_(good_component_locs, good_component_locs)]
                umi_support_subset = umis_per_seq[np.ix_(good_component_locs)]

                #calculate total number of mutations needed to explain data
                C = (D_subset*umi_support_subset).sum(axis=1)

                best_seq_idx = good_component_locs[np.argmin(C)]
                best_seq = unique_vdj_seqs[best_seq_idx]

                cell_vdj_calls.append(best_seq)

                for cindex in good_component_locs:
                    corrected_vdj_dict[it].update({unique_vdj_seqs[cindex] : best_seq})

            vdj_call_dict[it] = cell_vdj_calls

        else:
            is_cell_dict[it] = False
            is_multiplet_dict[it] = False
            vdj_call_dict[it] = None

    else:
        if max_umi_support > 1:
            is_cell_dict[it] = True
            is_multiplet_dict[it] = False
            vdj_call_dict[it] = [unique_vdj_seqs[-1]]
            corrected_vdj_dict[it][unique_vdj_seqs[-1]] = unique_vdj_seqs[-1]
        else:
            is_cell_dict[it] = False
            is_multiplet_dict[it] = False
            vdj_call_dict[it] = None


cb_df['is_cell'] = cb_df.index.map(lambda x: is_cell_dict[x]).astype(bool)
cb_df['is_multiplet'] = cb_df.index.map(lambda x: is_multiplet_dict[x]).astype(bool)
cb_df['vdj_calls'] = cb_df.index.map(lambda x: vdj_call_dict[x])
cb_df['multiplet'] = cb_df.vdj_calls.str.len()


print("total droplets per sample:", file=sys.stderr)
print(cb_df.groupby('sample_uid')['is_cell'].agg(sum), file=sys.stderr)

print("total cells per sample:", file=sys.stderr)
print(cb_df.groupby('sample_uid')['multiplet'].agg(sum), file=sys.stderr)

print("Multiplet distribution:", file=sys.stderr)

multiplet_distribution = cb_df.groupby(['sample_uid','multiplet']).size().reset_index()

multiplet_distribution = multiplet_distribution.pivot(index='multiplet',
                                                      columns='sample_uid',
                                                      values=0)

multiplet_distribution = multiplet_distribution/ multiplet_distribution.sum(axis=0)

print(multiplet_distribution, file=sys.stderr)

#collect info from ambient droplets
ambient_droplets = deepcopy(cb_df[(~cb_df['is_cell'].astype(bool))])
ambient_droplets['vdj_c_umi_r'] = ambient_droplets.apply(lambda x: [(x.vdj_sequence[i],
                                                                  x.c_gene_10X[i],
                                                                  str(x.umis[i]),
                                                                  str(x.reads[i])) for i in range(len(x.reads))
                                                                 ], axis=1)

ambient_droplets['droplet_id'] = ambient_droplets['sample_uid'] + "+" + ambient_droplets['cb'] + "+" + ambient_droplets['locus']
ambient_droplets = ambient_droplets[['droplet_id', 'vdj_c_umi_r']]
ambient_droplets = ambient_droplets.explode('vdj_c_umi_r')

called_cells = deepcopy(cb_df[cb_df['is_cell']])

ambient_vdj_dict = {}

cell_vdj_dict = {}
cell_c_call_dict = {}
cell_umi_count_dict = {}
cell_read_count_dict = {}

# for each called cell, assign isotype and record ambient RNA vdjs
UNIQUE_C_CALL = 0
MAJORITY_C_CALL_POSSIBLE = 0
AMBIGUOUS_C_CALL = 0

Ds = []
Cs = []

coencapsulated_ambient_rna = []

for it, row in called_cells.iterrows():
    # first record ambient RNA
    droplet_id = row.sample_uid + "+" + row.cb + "+" + row.locus

    ambient_vdj_info = [(row.vdj_sequence[it2],
                         row.c_gene_10X[it2],
                         row.umis[it2],
                         row.reads[it2])

                        for it2 in range(len(row.vdj_sequence))
                            if corrected_vdj_dict[it].get(row.vdj_sequence[it2],None) is None]
    droplet_ids = [droplet_id]*len(ambient_vdj_info)

    coencapsulated_ambient_rna.append(pd.DataFrame({'droplet_id':droplet_ids,
                                 'vdj_c_umi_r':ambient_vdj_info}))

    corrected_vdj_info = [(corrected_vdj_dict[it][row.vdj_sequence[it2]],
                           row.c_gene_10X[it2],
                           row.umis[it2],
                           row.reads[it2])

                        for it2 in range(len(row.vdj_sequence))
                            if not (corrected_vdj_dict[it].get(row.vdj_sequence[it2], None) is None)]

    # first collapse reads made collapsible by vdj error-correction
    cell_dict = {vdj:{} for vdj in row.vdj_calls}
    cell_c_call_dict[it] = {}
    cell_umi_count_dict[it] = {}
    cell_read_count_dict[it] = {}

    #find consensus C gene call
    for vdj, c_call, umis, reads in corrected_vdj_info:
        cell_dict[vdj][c_call] = (cell_dict[vdj].get(c_call, (0,0))[0] + int(umis),
                                  cell_dict[vdj].get(c_call, (0,0))[1] + int(reads))

    for vdj, umis_per_c_call in cell_dict.items():
        cell_umi_count_dict[it].update({vdj:{}})

        # first check whether c_call is unique
        if len(umis_per_c_call.keys()) == 1:
            UNIQUE_C_CALL += 1
            C_CALL = list(umis_per_c_call.keys())[0]
        else:
            #if not check whether best-supported C call
            #is supported by at least twice the number of UMIs as second-best C call

            total_umis = sum(u for u, r in umis_per_c_call.values())
            fractional_umis_per_c_call = sorted((umis_per_c_call[c][0]/total_umis, c)
                                            for c in umis_per_c_call.keys())

            best = fractional_umis_per_c_call[-1]
            second_best = fractional_umis_per_c_call[-2]

            if (second_best[0] / best[0] <= 0.5):
                MAJORITY_C_CALL_POSSIBLE += 1
                C_CALL = best[-1]

            else:
                AMBIGUOUS_C_CALL += 1
                C_CALL = 'ambiguous'


        cell_c_call_dict[it].update({vdj: C_CALL})

        #collect info about umi and read numbers for this cell
        cell_umi_count_dict[it][vdj] = sum(u for u, r in umis_per_c_call.values())
        cell_read_count_dict[it][vdj] = sum(r for u, r in umis_per_c_call.values())


print(f"Unique C gene call: {UNIQUE_C_CALL}", file=sys.stderr)
print(f"Majority C call: {MAJORITY_C_CALL_POSSIBLE}", file=sys.stderr)
print(f"Total confident C calls: {(UNIQUE_C_CALL + MAJORITY_C_CALL_POSSIBLE)}",
        file=sys.stderr)
print(f"Total cells called: {(UNIQUE_C_CALL + MAJORITY_C_CALL_POSSIBLE + AMBIGUOUS_C_CALL)}",
        file=sys.stderr)

called_cells['consensus_c_call'] = called_cells.index.map(cell_c_call_dict)
called_cells['n_umis'] = called_cells.index.map(cell_umi_count_dict)
called_cells['n_reads'] = called_cells.index.map(cell_read_count_dict)

called_cells['vdjc_info'] = called_cells.apply(lambda x: [(y, 
                                                           x.consensus_c_call[y], 
                                                           x.n_umis[y], 
                                                           x.n_reads[y])
                                                            for y in x.vdj_calls], axis=1)

called_cells = called_cells[['sample_uid', 'cb', 'locus', 'vdjc_info']]
called_cells = called_cells.explode('vdjc_info')
called_cells['vdj_call'] = called_cells.vdjc_info.map(lambda x: x[0])
called_cells['c_call'] = called_cells.vdjc_info.map(lambda x: x[1])
called_cells['n_umis'] = called_cells.vdjc_info.map(lambda x: x[2])
called_cells['n_reads'] = called_cells.vdjc_info.map(lambda x: x[3])

called_cells.to_csv(f'{OUTDIR}/{OUTNAME}_called_cells.tsv.gz', sep='\t', index=False)

print('Cell calling done!', file=sys.stderr)


print('Compiling information about ambient VDJ sequences...', file=sys.stderr)
coencapsulated_ambient_rna = pd.concat(coencapsulated_ambient_rna, ignore_index=True)
ambient_droplets = pd.concat([ambient_droplets, coencapsulated_ambient_rna], ignore_index=True)

ambient_droplets['vdj_sequence'] = ambient_droplets.vdj_c_umi_r.map(lambda x: x[0])
ambient_droplets['c_call'] = ambient_droplets.vdj_c_umi_r.map(lambda x: x[1])
ambient_droplets['umis'] = ambient_droplets.vdj_c_umi_r.map(lambda x: int(x[2]))
ambient_droplets['reads'] = ambient_droplets.vdj_c_umi_r.map(lambda x: int(x[3]))

ambient_droplets['sample_uid'] = ambient_droplets.droplet_id.map(lambda x: x.split("+")[0])
ambient_droplets['cb'] = ambient_droplets.droplet_id.map(lambda x: x.split("+")[1])
ambient_droplets['locus'] = ambient_droplets.droplet_id.map(lambda x: x.split("+")[2])

ambient_vdjs = ambient_droplets.groupby(['vdj_sequence', 'c_call', 'sample_uid', 'locus']).agg(list)
ambient_vdjs = ambient_vdjs.reset_index()

ambient_vdjs['n_droplets'] = ambient_vdjs.cb.map(lambda x: len(set(x)))
ambient_vdjs['n_umis'] = ambient_vdjs.umis.map(lambda x: sum(x))
ambient_vdjs['n_reads'] = ambient_vdjs.reads.map(lambda x: sum(x))

ambient_vdjs['cb'] = ambient_vdjs.cb.str.join(",")
ambient_vdjs['umis_per_droplet'] = ambient_vdjs.umis.map(lambda x: ",".join([str(y) for y in x]))
ambient_vdjs['reads_per_droplet'] = ambient_vdjs.reads.map(lambda x: ",".join([str(y) for y in x]))


ambient_vdjs.to_csv(f'{OUTDIR}/{OUTNAME}_ambient_vdjs.tsv.gz', sep='\t', index=False)
