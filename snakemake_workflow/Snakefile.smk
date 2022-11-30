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
# samplesheets = samplesheets[samplesheets.donor.str.contains('TBd3|TBd5|TBd6')]
samplesheets = samplesheets[samplesheets.sample_descriptor.str.contains('B200')]
# convenience variables
base = config["base"]
sample_uids = samplesheets.index.to_list()
donors = list(set(samplesheets[samplesheets.species == "human"].donor.to_list()))
os.makedirs(base, exist_ok=True)


rule all:
    input:
        "{}/analysis/scanpy/gex_object.h5ad".format(base),
    params:
        name="all",
        partition="quake",
    threads: 1


include: "rules/gex.smk"
include: "rules/vdjc.smk"


def samplesheet_lookup(idx, col):
    return samplesheets.loc[idx, col]


localrules:
    aggregate_h5ads,
