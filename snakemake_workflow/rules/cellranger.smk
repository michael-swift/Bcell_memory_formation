import sys
import os

# -----------------------------
# Configuration

species = config["species"]
### Cell Ranger

#### CellRanger VDJ #####
multi_config = "multi_config"


# construct multiconfig file from samplesheet
# done manually but with script rn
# rule cellranger mkfastq
rule cellranger_multi:
    input:
        workflow.basedir + "/config/multi_config_{sample_uid}.csv",
    output:
        config="{base}/10X/{sample_uid}/outs/per_sample_outs/{sample_uid}/web_summary.html",
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"]
    resources: mem_mb="210000",partition="quake",disk_mb="8000",time="6-0"
    threads:20
    shell:
        "mkdir -p {params.base}/10X &&"
        " cd {params.base}/10X &&"
        " rm -rf {wildcards.sample_uid} && {params.cell_ranger}/cellranger multi --id={wildcards.sample_uid} --csv={input}"
