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

samplesheets["sample_uid"] = (
    samplesheets["donor"]
    + "_"
    + samplesheets["sample_type"]
    + "_"
    + samplesheets["sample_descriptor"]
)
samplesheets["species"] = "human"
samplesheets.set_index("sample_uid", inplace=True)
# convenience variables
base = config["base"]
donors = list(set(samplesheets[samplesheets.species == "human"].donor.to_list()))
os.makedirs(base, exist_ok=True)

samplesheets['expected_cells'] = samplesheets['expected_cells_thousands'].astype(int) * 1000


#samplesheets = samplesheets[samplesheets.donor != "TBd5"]
sample_uids = samplesheets.index.to_list()
rule all:
    input:
        expand("{base}/per_sample/cellranger_vdj/{sample_uid}/outs/web_summary.html", base = base, sample_uid = sample_uids),
        "{}/analysis/scanpy/gex_object.h5ad.gz".format(base),
    params:
        name="all",
        partition="quake",
    threads: 1


include: "rules/gex.smk"
include: "rules/vdjc.smk"
#include: "rules/cellranger.smk"

def samplesheet_lookup(idx, col):
    return samplesheets.loc[idx, col]

#localrules:preprocess_scanpy
