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
samplesheets["expected_cells"] = (
    samplesheets["expected_cells_thousands"].astype(int) * 1000
)
samplesheets["species"] = "human"
donors = list(set(samplesheets[samplesheets.species == "human"].donor.to_list()))

# filter samplesheet:
samplesheets = samplesheets[samplesheets.sample_uid.str.contains('TBd4_frozen_PB')]
# samplesheets = samplesheets.iloc[:2]
samplesheets_gex = samplesheets[samplesheets.lib_type == "gex"]
samplesheets_gex.set_index("sample_uid", inplace=True)
samplesheets_vdj = samplesheets[samplesheets.lib_type == "vdj"]
# set wildcards
sample_uids_vdj = samplesheets_vdj.sample_uid.to_list()
samplesheets = samplesheets_gex
sample_uids = samplesheets.index.to_list()
print(sample_uids)
# make processsed data dir
os.makedirs(base, exist_ok=True)


rule all:
    input:
        #expand("{base}/per_sample/cellranger_vdj/{sample_uid_vdj}/outs/web_summary.html", base = base['gex'], sample_uid_vdj = sample_uids_vdj),
        "{}/analysis/scanpy/gex_object.h5ad.gz".format(base['gex']),
        #expand("{base}/per_sample/fastqc/{sample_uid_vdj}/", base = base['gex'], sample_uid_vdj = sample_uids_vdj),
        #expand("{base}/per_sample/star_solo_vdj/{sample_uid_vdj}/Aligned.out.bam", base = base['gex'], sample_uid_vdj = sample_uids_vdj),
        expand("{base}/aggregated/lineage_clustering/final_lineage_ids/{donor}.tsv.gz", base=base['vdj'], donor=donors),
        expand("{base}/aggregated/vtrees/{which}/{donor}_v_trees.tsv", base = base['vdj'], which = ['cells'], donor = donors),
        expand("{base}/aggregated/cell_calls/{donor}_called_cells_vdj_annotated.tsv.gz", base = base['vdj'], donor = donors)
    params:
        name="all",
        partition="quake",
    threads: 1


include: "rules/gex.smk"
include: "rules/vdjc.smk"
include: "rules/qc.smk"
include: "rules/get_resources.smk"
include: "rules/cell_calling.smk"


def samplesheet_lookup(idx, col):
    return samplesheets.loc[idx, col]
