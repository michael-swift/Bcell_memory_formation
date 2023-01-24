import os
import pandas as pd
import glob
from collections import defaultdict

shell.prefix("set +euo pipefail;")
configfile: "config/path_config.yaml"
base = config["base"]
library_info = config["library_info"]
libuid_to_datadirs = pd.read_table(library_info, sep="    ", engine="python")
libuid_to_datadirs = libuid_to_datadirs.set_index("library_uid")["path"].to_dict()

# TODO put in config
samplesheets = pd.concat(
    [
        pd.read_table("{}".format(x), dtype="str", sep="\t", engine="python")
        for x in config["samplesheets"]
    ],
    ignore_index=True,
)

# Parse Samplesheet
samplesheets['expected_cells'] = samplesheets['expected_cells_thousands'].astype(float) * 1000
samplesheets["species"] = "human"
donors = list(set(samplesheets[samplesheets.species == "human"].donor.to_list()))


samplesheets_gex = samplesheets[samplesheets.lib_type == 'gex']
samplesheets_gex.set_index("sample_uid", inplace=True)
samplesheets_vdj = samplesheets[samplesheets.lib_type == 'vdj']
# set wildcards
sample_uids_vdj = samplesheets_vdj.sample_uid.to_list()
samplesheets = samplesheets_gex
sample_uids = samplesheets.index.to_list()
# make processsed data dir
os.makedirs(base, exist_ok=True)


rule all:
    input:
        expand("{base}/aggregated/lineage_clustering/final_lineage_ids/{donor}.tsv.gz", base=base, donor=donors),
        expand("{base}/aggregated/vtrees/pseudobulk/{donor}_v_trees.tsv", base = base, donor = donors),
        expand("{base}/aggregated/cell_calls/{donor}_called_cells_vdj_annotated.tsv.gz", base = base, donor = donors)
        #"{}/analysis/scanpy/gex_object.h5ad.gz".format(base),
    params:
        name="all",
        partition="quake",
    threads: 1


include: "rules/gex.smk"
include: "rules/vdjc.smk"
include: "rules/cell_calling.smk"

def samplesheet_lookup(idx, col):
    return samplesheets.loc[idx, col]

#localrules:preprocess_scanpy
