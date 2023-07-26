import os
import pandas as pd
import glob
from collections import defaultdict

shell.prefix("set +euo pipefail;")
configfile: "config/path_config.yaml"
base = config["base"]
# TODO put in config
samplesheets = pd.concat([pd.read_table("{}".format(x), dtype="str", sep="\t", engine="python") for x in config["samplesheets"]], ignore_index=False)
# Parse Samplesheet
samplesheets["expected_cells"] = (
    samplesheets["expected_cells_thousands"].astype(int) * 1000
)

samplesheets["species"] = "human"
donors = list(set(samplesheets[samplesheets.species == "human"].donor.to_list()))
samplesheets_gex = samplesheets[samplesheets.lib_type == "gex"]
samplesheets_gex.set_index("sample_uid", inplace=True)
samplesheets_vdj = samplesheets[samplesheets.lib_type == "vdj"]
# setup wildcards
sample_uids_vdj = samplesheets_vdj.sample_uid.to_list()
samplesheets = samplesheets_gex
sample_uids = samplesheets.index.to_list()
tissues = samplesheets.tissue.to_list()
print(sample_uids)

# make processsed data dir
for k in base.keys():
    os.makedirs(base[k], exist_ok=True)


rule all:
    input:
        #expand("{base}/per_sample/cellranger_vdj/{sample_uid_vdj}/outs/web_summary.html", base = base['gex'], sample_uid_vdj = sample_uids_vdj),
        #expand("{base_gex}/annotate/gex_object.h5ad.gz", base_gex = base['gex']),
        #expand("{base}/per_sample/fastqc/{sample_uid_vdj}/", base = base['gex'], sample_uid_vdj = sample_uids_vdj),
        #expand("{base}/per_sample/star_solo_vdj/{sample_uid_vdj}/Aligned.out.bam", base = base['gex'], sample_uid_vdj = sample_uids_vdj),
        expand(
            "{base_gex}/outs/{celltypes}.h5ad.gz",
            base_gex=base["gex"],
            celltypes=["ASC", "MB", "NB_other", "bcells", "only_igh", "all_cells"],
        ),
        #expand("{base}/aggregated/vdj/{donor}_combined.tsv.gz", base=base['vdj'], donor=donors),
        #expand("{base}/aggregated/lineage_clustering/final_lineage_ids/{donor}.tsv.gz", base=base['vdj'], donor=donors),
        #expand("{base}/aggregated/vtrees/{which}/{donor}_v_trees.tsv", base = base['vdj'], which = ['cells'], donor = donors),
        expand("{base}/all_vdj_cell_calls_IGH.tsv.gz", base = base['vdj']),
        expand("{base}/integrated_cell_calls.tsv.gz", base=base['vdj'])
    params:
        name="all",
        partition="quake",
    resources:
        threads=1,


include: "rules/gex.smk"
include: "rules/get_resources.smk"
include: "rules/annotate.smk"
include: "rules/qc.smk"
include: "rules/cell_calling.smk"
include: "rules/associate_gex_vdj.smk"

# localrules:merge_vdj
def samplesheet_lookup(idx, col):
    return samplesheets.loc[idx, col]
