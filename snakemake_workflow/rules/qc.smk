from glob import glob


def get_r1_fqgz_using_wildcards(wildcards):
    """ Returns R1 and R2 gzipped fastq file paths as a list """

    r1 = glob('{fastqs}/{sample_uid_vdj}_*R1*.fastq.gz'.format(fastqs=config['vdj_fastq_dir'],
                                                     sample_uid_vdj=wildcards['sample_uid_vdj']
                                                     ))

    r2 = glob('{fastqs}/{sample_uid_vdj}_*R2*.fastq.gz'.format(fastqs=config['vdj_fastq_dir'],                                                   sample_uid_vdj=wildcards['sample_uid_vdj']                                                        ))
    return r1

def get_r2_fqgz_using_wildcards(wildcards):
    """ Returns R1 and R2 gzipped fastq file paths as a list """

    r2 = glob('{fastqs}/{sample_uid_vdj}_*R2*.fastq.gz'.format(fastqs=config['vdj_fastq_dir'],                                                   sample_uid_vdj=wildcards['sample_uid_vdj']                                                        ))
    return r2

rule fastqc:
    output:
        touch("{base}/per_sample/fastqc/{sample_uid_vdj}/done.txt"),
        directory("{base}/per_sample/fastqc/{sample_uid_vdj}/")
    log:
        "{base}/logs/{sample_uid_vdj}/fastqc.log",
    params:
        name="fastq_qc",
        base=config["base"],
        fastq_dir=config["vdj_fastq_dir"],
    resources:
        mem_mb="32000",
        partition="quake,owners,normal",
        disk_mb="8000",
        time="0-4",
        threads=4
    threads: 4
    conda: config["workflow_dir"] + "/envs/single-cell.yaml"
    shell:
        "mkdir -p {output[1]} && fastqc -t 8 {params.fastq_dir}{wildcards.sample_uid_vdj}_* -o {output[1]}"
