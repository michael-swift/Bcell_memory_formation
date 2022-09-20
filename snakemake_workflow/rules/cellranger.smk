import sys
import os

# -----------------------------
# Configuration

species = config["species"]
### Cell Ranger

#### CellRanger VDJ #####
multi_config = "multi_config"

rule cellranger_count:
    input:
    output:
        config="{base}/10X/{sample_uid}/outs/web_summary.html",
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="count_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        transcriptome=config["transcriptome"],
        fastq_dir=config['fastq_dir']
    resources: mem_mb="210000",partition="quake",disk_mb="8000",time="2-0"
    threads:20
    shell:
        "mkdir -p {params.base}/10X &&"
        " cd {params.base}/10X &&"
        " rm -rf {wildcards.sample_uid} && {params.cell_ranger}/cellranger count --id={wildcards.sample_uid} --transcriptome={params.transcriptome} --fastqs {params.fastq_dir} --sample={wildcards.sample_uid} --localcores=20 --localmem=210 --nosecondary=True --no-bam=True"
