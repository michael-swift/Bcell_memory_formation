rule fastqc:
    output:
        "{base}/per_sample/fastqc/{sample_uid_vdj}/done.txt",
        directory("{base}/per_sample/fastqc/{sample_uid_vdj}/")
    log:
        "{base}/logs/{sample_uid_vdj}/fastqc.log",
    params:
        name="fastq_qc",
        base=config["base"],
        fastq_dir=config["vdj_fastq_dir"],
    resources:
        mem_mb="32000",
        partition="quake,owners",
        disk_mb="8000",
        time="0-3",
    threads: 20
    conda: config["workflow_dir"] + "/envs/single-cell.yaml"
    shell:
        "mkdir -p {output[1]} && fastqc {params.fastq_dir}{wildcards.sample_uid_vdj}_* -o {output[1]} && touch {output[0]}"
