rule cellranger_vdj:
    output:
        "{base}/per_sample/cellranger_vdj/{sample_uid}/outs/web_summary.html",
    log:
        "{base}/logs/{sample_uid}/cellranger.log",
    params:
        name="vdj_cellranger",
        base=config["base"],
        cell_ranger=config["cell_ranger"],
        vdj_reference=config["vdj_reference"],
        fastq_dir=config["vdj_fastq_dir"],
        inner_primers=config["inner_primers"]
    resources:
        mem_mb="120000",
        partition="quake,owners,normal",
        disk_mb="8000",
        time="1-0",
    threads: 20
    shell:
        "mkdir -p {params.base}/cellranger_vdj && "
        "cd {params.base}/cellranger_vdj && "
        "rm -rf {wildcards.sample_uid} && "
        "{params.cell_ranger}/cellranger vdj --id={wildcards.sample_uid} "
        "--reference={params.vdj_reference} --fastqs {params.fastq_dir} "
        "--sample={wildcards.sample_uid} --localcores=20 > {log}"
