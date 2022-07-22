rule dandelion_image:
    """ pull the immcantation image using singularity

    """
    output:
        "{}/resources/repertoire_analysis/{}.sif".format(workflow.basedir, "sc-dandelion"),
    threads: 1
    params:
        name="container_pull",
        docker_address="docker://immcantation/suite:4.1.0",
    resources:
        mem_mb=10000,
    shell:
        "module load system && "
        "mkdir -p $(dirname {output}) && cd $(dirname {output}) && "
        "img=$(basename {output}) && "
        "singularity pull library://kt16/default/sc-dandelion:latest"


rule get_transcriptome_ref:
    output:
        transcriptome="transcriptome".format(workflow.basedir),
        outdir=directory("{}/db/transcriptome/".format(workflow.basedir)),
    params:
        ref="{}/db/10X/{}/".format(workflow.basedir, config["vdj_ref_prefix"]),
        outdir="{}/db/10X/".format(workflow.basedir),
        species="Homo sapiens",
        link=config["10X_transcriptome"],
        name=config["10X_transcriptome"].split("/")[-1]
    conda:
        "../envs/cellranger.yaml"
    shell:
        "mkdir -p {output.outdir} "
        "&& cd {output.outdir} "
        "&& wget {params.link} "
        "&& tar -xvf {params.name}"
