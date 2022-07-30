from pacbio_vdj_utils.trie_dict_utils import *

import numpy as np
import pandas as pd

import argparse
import sys, time
#####################################################    

parser = argparse.ArgumentParser(description='Whitelists cell barcodes based on provided whitelist '
                                            'and returns a list of reads and UMIs associated with '
                                            'whitelisted barcodes')
parser.add_argument('cb_umi_seqid_map')
parser.add_argument('-whitelist', required=True, 
                    help='path to cell barcode whitelist')
parser.add_argument('-outdir', default='.', help="output directory (defualt: working directory)")
args = parser.parse_args()

CB_UMI_SEQ_ID_MAP_FILE = args.cb_umi_seqid_map
WHITELIST_FILE = args.whitelist
outdir = args.outdir

#####################################################

WHITELIST = pd.read_table(WHITELIST_FILE, header=None)[0].values

start = time.time()
trie = read_trie(WHITELIST)
end = time.time()
sys.stderr.write("INFO: Constructing whitelist trie took %gs\n" % (end - start))

#####################################################

seq_id_df = pd.read_table(CB_UMI_SEQ_ID_MAP_FILE)

seq_id_df['CBs'] = seq_id_df.SEQ_ID_INFO.map(lambda x: [y for y in x.split("|")[0].split("=")[-1].split(",")])
seq_id_df['CBs'] = seq_id_df.CBs.map(lambda x: [y.split("_")[0] for y in x])
all_cell_barcodes = seq_id_df['CBs'].explode().unique()

CBs = {k:None for k in all_cell_barcodes}

start = time.time()

exact = 0
corrected = 0
total = 0

for cb in CBs.keys():
    total += 1
    if len(search(cb, trie, 0)) == 1:
        CBs[cb] = cb
        exact += 1
    else:
        single_edit = search(cb, trie, 1)
        if len(single_edit) == 1:
            if single_edit[0][0] in CBs.keys():
                CBs[cb] = cb
                corrected += 1

end = time.time()

sys.stderr.write("INFO: Whitelisting and error correction took %gs\n" % (end - start))

sys.stderr.write("INFO: {} barcodes processed\n".format(total))
sys.stderr.write("INFO: {} barcodes whitelisted ({} exact)\n".format((exact + corrected),exact))
sys.stderr.write("INFO: {} cell barcodes discarded\n".format((total - exact - corrected)))


sample_uid=CB_UMI_SEQ_ID_MAP_FILE.split("/")[-1].split("_")[0]

outwriter = open("{}/{}_consensus-pass_seq_ids_barcodes_whitelisted.tsv".format(outdir,sample_uid), 'w')
outwriter.write("\t".join(['SEQ_ID', 'SEQ_ID_INFO']) + "\n")

start = time.time()

discarded_seq_ids = 0

for it, row in seq_id_df.iterrows():
    CBs_UMIs = row.SEQ_ID_INFO.split("|")[0].split("=")[-1].split(",")

    cbs = [x.split("_")[0] for x in CBs_UMIs]
    umis = [x.split("_")[1] for x in CBs_UMIs]
    reads = row.SEQ_ID_INFO.split("|")[1].split("=")[-1].split(",")

    NEW_SEQ_ID_INFO = "BARCODE="
    NEW_SEQ_ID_INFO += ",".join(["_".join([CBs[cbs[it]], umis[it]]) 
                                for it in range(len(cbs)) if not(CBs[cbs[it]] is None)])
    NEW_SEQ_ID_INFO += "|CONSCOUNT="
    NEW_SEQ_ID_INFO += ",".join([reads[it]
                                for it in range(len(cbs)) if not(CBs[cbs[it]] is None)])

    if NEW_SEQ_ID_INFO.startswith("BARCODE=|CONSCOUNT"):
        discarded_seq_ids += 1
    else:
        outwriter.write("\t".join([str(row.SEQ_ID), NEW_SEQ_ID_INFO]) + "\n")

outwriter.close()
sys.stderr.write("INFO: Writing new seq_id map took {}s\n".format((time.time()-start)))

sys.stderr.write("INFO: {} seqids have been discarded due to no valid CBs\n".format(discarded_seq_ids))

