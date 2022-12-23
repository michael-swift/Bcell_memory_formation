rule cellranger_vdj:
    output:
        "{base}/per_sample/cellranger_vdj/{sample_uid_vdj}/outs/web_summary.html",
    log:
        "{base}/logs/{sample_uid_vdj}/cellranger.log",
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
        time="0-12",
    threads: 20
    shell:
        "mkdir -p {params.base}/per_sample/cellranger_vdj && "
        "cd {params.base}/per_sample/cellranger_vdj && "
        "rm -rf {wildcards.sample_uid_vdj} && "
        "{params.cell_ranger}/cellranger vdj --id={wildcards.sample_uid_vdj} "
        "--reference={params.vdj_reference} --fastqs {params.fastq_dir} "
        "--sample={wildcards.sample_uid_vdj} --localcores=20 > {log}"
