import os
import pandas as pd
import glob
from collections import defaultdict

shell.prefix("set +euo pipefail;")
configfile: "config/config.yaml"
base = config["base"]
library_info = config["library_info"]
libuid_to_datadirs = pd.read_table(library_info, sep="    ", engine="python")
libuid_to_datadirs = libuid_to_datadirs.set_index("library_uid")["path"].to_dict()

# TODO put in config
samplesheets = pd.concat(
    [
        pd.read_table("{}".format(x), dtype="str", sep="    ", engine="python")
        for x in config["samplesheets"]
    ],
    ignore_index=True,
)

# print(samplesheets.columns)
samplesheets["sample_uid"] = samplesheets["sample_id"]
samplesheets.set_index("sample_uid", inplace=True)

# convenience variables
base = config["outs_basedir"]
sample_uids = samplesheets.index.to_list()
donors = list(set(samplesheets[samplesheets.species == "human"].donor.to_list()))
os.makedirs(base, exist_ok=True)


rule all:
    input:
        expand(
            "{base}/10X/{sample_uid}/outs/per_sample_outs/{sample_uid}/web_summary.html",
            base=config["base"],
            sample_uid=sample_uids,
        ),
    params:
        name="all",
        partition="quake",
    threads: 1


include: "rules/get_resources.smk"
include: "rules/cellranger.smk"
include: "rules/vdjc.smk"


def samplesheet_lookup(idx, col):
    return samplesheets.loc[idx, col]
