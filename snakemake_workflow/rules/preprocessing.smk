def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    sample_uid = wildcards.sample_uid
    # get sample_index based and sample_uid
    sample_index = samplesheets.loc[sample_uid]["sample_index"]
    # map libuid to actual sherlock path
    libuid = int(samplesheets.loc[sample_uid, "library_uid"])
    data_path = libuid_to_datadirs[libuid]
    fastq = (
        data_path
        + "/"
        + "demultiplex.{}--{}.hifi_reads.fastq.gz".format(sample_index, sample_index)
    )
    return fastq


def get_lib_type(wildcards):
    sample_uid = wildcards.sample_uid
    return samplesheets.loc[sample_uid, "lib_type"]


rule filter_quality:
    input:
        get_fastq,
    output:
        temp_fq=temp("{base}/filter_quality/{sample_uid}.fastq"),
        fastq="{base}/filter_quality/{sample_uid}_quality-pass.fastq.gz",
        logfile="{base}/prestologs/{sample_uid}_filter-quality.log",
    log:
        "{base}/logs/{sample_uid}_filter-quality.log",
    params:
        output_dir=config["base"],
    threads: 16
    conda:
        "../envs/presto.yaml"
    shell:
        """
        zcat {input} > {output.temp_fq}

        FilterSeq.py quality \
        -s {output.temp_fq} \
        -q 20 \
        --outdir {params.output_dir}/filter_quality \
        --outname {wildcards.sample_uid} \
        --log {output.logfile} \
    --nproc {threads} \
        > {log}

        gzip {base}/filter_quality/{wildcards.sample_uid}_quality-pass.fastq
        """


rule fix_read_orientation:
    input:
        "{base}/filter_quality/{sample_uid}_quality-pass.fastq.gz",
    output:
        fastq="{base}/fix_read_orientation/{sample_uid}_quality-pass_cleaned.fastq.gz",
    log:
        "{base}/logs/{sample_uid}_flip.log",
    params:
        scripts=config["scripts"],
    conda:
        "../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/fix_read_orientation.py {input} "
        "-outdir {wildcards.base}/fix_read_orientation "
        "2> {log}"


rule parse_reads:
    input:
        "{base}/fix_read_orientation/{sample_uid}_quality-pass_cleaned.fastq.gz",
    output:
        amplicon="{base}/parse_reads/{sample_uid}_quality-pass_amplicon.fastq.gz",
        bc="{base}/parse_reads/{sample_uid}_quality-pass_bc.fastq.gz",
    log:
        "{base}/logs/{sample_uid}_parse.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
        read_config=config["read_config"],
        lib_type=get_lib_type,
        minibulk=lambda wildcards: "--minibulk "
        if "minibulk" in wildcards.sample_uid
        else "",
    shell:
        "python {params.scripts}/parse_reads_stringent_10X.py {input} "
        "-config {params.read_config} "
        "-outdir {wildcards.base}/parse_reads "
        "--{params.lib_type} "
        "{params.minibulk}"
        "2> {log}"


rule call_cb_whitelist:
    input:
        rules.parse_reads.output.bc,
    output:
        "{base}/whitelists/{sample_uid}_cb_whitelist.txt",
    conda:
        "../envs/umi_tools.yaml"
    log:
        "{base}/logs/{sample_uid}_cb_whitelist.log",
    shell:
        "umi_tools whitelist -I {input} "
        "--bc-pattern='(?P<cell_1>.{{16}})(?P<umi_1>.{{10}})' "
        "--knee-method=density "
        "--extract-method=regex "
        "--plot-prefix={wildcards.base}/plots/{wildcards.sample_uid} "
        "-L {log} "
        "> {output} "


rule align_umi_groups:
    input:
        fastq=rules.parse_reads.output.amplicon,
    output:
        tmp_fq=temp("{base}/align_umi_groups/{sample_uid}.fastq"),
        fastq="{base}/align_umi_groups/{sample_uid}_align-pass.fastq",
        logfile="{base}/prestologs/{sample_uid}_align-sets.log",
    log:
        "{base}/logs/{sample_uid}_align.log",
    threads: 20
    conda:
        "../envs/presto.yaml"
    params:
        scripts=config["scripts"],
    shell:
        """
        zcat {input} > {output.tmp_fq}

        AlignSets.py muscle \
            -s {output.tmp_fq} \
            --outname {wildcards.sample_uid} \
            --outdir {wildcards.base}/align_umi_groups \
            --log {output.logfile} \
            > {log}
        """


rule build_consensus:
    input:
        "{base}/align_umi_groups/{sample_uid}_align-pass.fastq",
    output:
        fastq="{base}/build_consensus/{sample_uid}_consensus-pass.fastq",
    log:
        prestolog="{base}/prestologs/{sample_uid}_build-consensus.log",
        runlog="{base}/logs/{sample_uid}_build-consensus.log",
    threads: 8
    conda:
        "../envs/presto.yaml"
    shell:
        "BuildConsensus.py "
        "-s {input} "
        "--bf BARCODE "
        "--cf BARCODE CPRIMER "
        "--act majority set "
        "--pf CPRIMER "
        "--prcons 0.6 "
        "--maxerror 0.1 "
        "--maxgap 0.5 "
        "--outname {wildcards.sample_uid} "
        "--outdir {wildcards.base}/build_consensus/ "
        "--log {log.prestolog} "
        "> {log.runlog}"


rule deduplicate_sample:
    input:
        "{base}/build_consensus/{sample_uid}_consensus-pass.fastq",
    output:
        "{base}/deduplicate_sample/{sample_uid}_consensus-pass_deduplicated.fasta",
        "{base}/deduplicate_sample/{sample_uid}_consensus-pass_seq_ids_barcodes.tsv",
    log:
        "{base}/logs/{sample_uid}_deduplicate.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/deduplicate_sequences_for_mapping.py "
        "{input} "
        "-outdir {wildcards.base}/deduplicate_sample "
        "-copy_fields BARCODE CONSCOUNT "
        "2> {log}"
