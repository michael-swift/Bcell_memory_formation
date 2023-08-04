smplesheet = samplesheets
include: "gex.smk"
include: "annotate.smk"
include: "cell_calling.smk"

rule integrate_gex_and_vdj_data:
    input:
        vdj=expand("{base_vdj}/all_vdj_cell_calls_IGH.tsv.gz", base_vdj = base['vdj']),
        gex=expand("{base_gex}/annotate/all_adata.obs.tsv.gz", base_gex = base['gex'])
    output:
        expand("{base_vdj}/integrated_cell_calls.tsv.gz", base_vdj = base['vdj'])
    params:
        scripts=config["vdj_scripts"],
        info=config['sample_relationships'],
    log:
        expand("{base_vdj}/logs/integrate_gex_and_vdj_data.log", base_vdj = base['vdj']),
    resources:
        mem_mb="128000",
    conda:
        "../envs/pacbio.yaml",
    shell:
        "python {params.scripts}/integrate_gex_and_vdj.py "
        "-vdj {input.vdj} "
        "-gex {input.gex} "
        "-output {output} "
        "-info {params.info} "
        "2> {log}"

rule annotate_gex_ambient:
    input:
        "{base_vdj}/integrated_cell_calls.tsv.gz",
    output:
        "{base_vdj}/integrated_cell_calls_ambient_annotated.tsv.gz",
    params:
        scripts=config["vdj_scripts"],
    log:
        "{base_vdj}/logs/annotate_ambient.log",
    resources:
        mem_mb="128000",
    conda:
        "../envs/pacbio.yaml",
    shell:
        "python {params.scripts}/annotate_ambient_rna_in_GEX_samples.py "
        "-input {input} "
        "-output {output} "
        "2> {log}"

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint prepare_for_olga:
    input:
        "{base_vdj}/integrated_cell_calls.tsv.gz",
    output:
        groups=directory("{base_vdj}/aggregated/olga/cdr3"),
        complete="{base_vdj}/aggregated/olga/cdr3/cdr3_compile.complete",
    log:
        "{base_vdj}/logs/prepare_for_olga.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
        base=base['vdj'],
    resources:
        mem_mb="64000",
        time="4:00:00",
    shell:
        """
        mkdir -p {params.base}/aggregated/olga/cdr3
        mkdir -p {params.base}/aggregated/olga/pgen
        
        python {params.scripts}/prepare_for_olga.py \
        -input {input} \
        -outdir {params.base}/aggregated/olga/cdr3 \
        2> {log}
        
        touch {output.complete}
        """

rule run_olga:
    input:
        "{base_vdj}/aggregated/olga/cdr3/cdr3nt_chunk-{group}.tsv",
    output:
        "{base_vdj}/aggregated/olga/pgen/pgen_marginal-{group}.tsv",
    log:
        "{base_vdj}/logs/olga/{group}.log",
    conda:
        "../envs/olga.yaml"
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="8000",
        time="8:00:00",
    shell:
        "olga-compute_pgen --humanIGH "
        "-i {input} "
        "-o {output} "
        "--seq_in 1 "
        "--display_off"

def aggregate_olga_output(wildcards):
    checkpoint_output = checkpoints.prepare_for_olga.get(
        **wildcards
    ).output.groups
    
    files = expand(
        "{base_vdj}/aggregated/olga/pgen/pgen_marginal-{group}.tsv",
        base_vdj=wildcards.base_vdj,
        group=glob_wildcards(
            os.path.join(checkpoint_output, "cdr3nt_chunk-{group}.tsv")
            ).group,
        )
    return files

rule concatenate_olga_output:
    input:
        files=aggregate_olga_output,
        checkpoint_complete="{base_vdj}/aggregated/olga/cdr3/cdr3_compile.complete",
    output:
        "{base_vdj}/pgen_marginal.tsv.gz",
    log:
        "{base_vdj}/logs/olga/concatenate_olga_output.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="131000",
        time="12:00:00",
    shell:
        "python {params.scripts}/concatenate_olga_output.py "
        "-input {input.files} "
        "-output {output} "
        "2> {log}"



