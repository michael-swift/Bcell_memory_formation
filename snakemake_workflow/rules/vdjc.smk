rule cellranger_vdj:
    output:
        csv="{base}/per_sample/cellranger_vdj/{sample_uid_vdj}/outs/all_contig_annotations.csv",
        fasta="{base}/per_sample/cellranger_vdj/{sample_uid_vdj}/outs/all_contig.fasta",
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
        partition="quake,owners",
        disk_mb="8000",
        time="0-16",
    threads: 20
    shell:
        "mkdir -p {params.base}/per_sample/cellranger_vdj && "
        "cd {params.base}/per_sample/cellranger_vdj && "
        "rm -rf {wildcards.sample_uid_vdj} && "
        "{params.cell_ranger}/cellranger vdj --id={wildcards.sample_uid_vdj} "
        "--reference={params.vdj_reference} --fastqs {params.fastq_dir} "
        "--sample={wildcards.sample_uid_vdj} --localcores=20 > {log}"

rule igblast:
    input:
        rules.cellranger_vdj.output.fasta,
    output:
        "{base}/per_sample/vdj_preprocess/{sample_uid_vdj}/igblast.tsv",
    conda:
        "../envs/pacbio.yaml",
    threads:
        16,
    params:
        organism="human",
        IGDBDIR=config["IGDBDIR"],
    shell:
        """
        ml load system libuv

        export IGDATA={params.IGDBDIR}
        
        igblastn \
                -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_ig_v \
                -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_ig_d \
                -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_ig_j \
                -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
                -domain_system imgt \
                -ig_seqtype Ig \
                -organism {params.organism} \
                -outfmt 19 \
                -query {input} \
                -out {output} \
                -num_threads {threads}
        """

rule filter_and_annotate:
    input:
        tenX=rules.cellranger_vdj.output.csv,
        igblast=rules.igblast.output,
    output:
        "{base}/per_sample/vdj_preprocess/{sample_uid_vdj}/igblast_filtered_annotated.tsv.gz",
    conda:
        "../envs/pacbio.yaml",
    log:
        "{base}/logs/{sample_uid_vdj}_filter_vdj.log",
    resources:
        mem_mb="131000",
    params:
        scripts=config["vdj_scripts"],
    shell:
        "python {params.scripts}/filter_igblast_and_annotate.py "
        "{input.igblast} "
        "-TXG {input.tenX} "
        "-outdir {wildcards.base}/per_sample/vdj_preprocess/{wildcards.sample_uid_vdj}  "
        "--verbose "
        "> {log} "
        "2> {log}"

def aggregate_input(wildcards):
    samplelist = ["{}/per_sample/vdj_preprocess/{}/igblast_filtered_annotated.tsv.gz".format(wildcards.base, sample_uid_vdj) for sample_uid_vdj in samplesheets_vdj[(samplesheets_vdj.donor == wildcards.donor)].sample_uid]
    return samplelist

def get_cb_umi_maps(wildcards):
    maplist = ["{}/whitelisted/{}_seq_ids_barcodes_whitelisted.tsv".format( config["outs_basedir"], sid) for sid in samplesheets[(samplesheets.donor==wildcards.donor) & (samplesheets.primer_set == wildcards.primer_set)].index]
    return maplist


rule combine_samples:
    input:
        aggregate_input,
    output:
        "{base}/aggregated/vdj/{donor}_combined.tsv.gz",
    conda:
        "../envs/pacbio.yaml",
    log:
        "{base}/logs/{donor}_combined.log",
    resources:
        mem_mb="131000",
    params:
        scripts=config["vdj_scripts"],
        samplesheet=config["samplesheet_vdj"],
    shell:
        "python {params.scripts}/aggregate_samples_by_donor.py "
        "{input} "
        "-samplesheet {params.samplesheet} "
        "-outdir {wildcards.base}/aggregated/vdj "
        "-name {wildcards.donor} "

rule call_germlines:
    input:
        rules.combine_samples.output,
    output:
        germlines="{base}/aggregated/grmlin/{donor}_combined_v_germlines.tsv",
        preprocessed="{base}/aggregated/grmlin/{donor}_combined_preprocessed.tsv.gz",
    log:
        "{base}/logs/grmlin/{donor}_grmlin.log",
    conda:
        "../envs/grmlin.yaml"
    params:
        organism=lambda wildcards: samplesheets[
                samplesheets.donor == str(wildcards.donor)]["species"].values[0],
        grmlin=config["grmlin"],
        IGDBDIR=config["IGDBDIR"],
    resources:
        mem_mb="131000",
        time="24:00:00",
    shell:
        "{params.grmlin}/grmlin "
        "{input} "
        "-annotate {params.IGDBDIR}/database/imgt_{params.organism}_ig_v "
        "-outdir {wildcards.base}/aggregated/grmlin "
        "--verbose "
        "-max_sequences 500000 "
        "> {log}"

# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint prepare_cdr3_groups_for_distance_evaluation:
    input:
        rules.call_germlines.output.preprocessed,
    output:
        groups=directory("{base}/aggregated/lineage_clustering/cdr3/{donor}/"),
    wildcard_constraints:
        donor="TBd[1-6]",
    log:
        "{base}/logs/cluster_lineages/{donor}_prepare_for_clustering.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="64000",
        time="24:00:00",
    shell:
        """
        mkdir -p {wildcards.base}/aggregated/lineage_clustering/cdr3/{wildcards.donor} 
        
        python {params.scripts}/prepare_distance_matrices.py \
        {input} \
        -outdir {wildcards.base}/aggregated/lineage_clustering/cdr3/{wildcards.donor} \
        -samplename {wildcards.donor} \
        2> {log}
        
        """

rule fasta_to_hamming_distance_matrix:
    input:
        "{base}/aggregated/lineage_clustering/cdr3/{donor}/{group}.fasta",
    output:
        "{base}/aggregated/lineage_clustering/cdr3/{donor}/{group}.npy",
    log:
        "{base}/logs/cluster_lineages/matrix_calc/cdr3/{donor}_{group}.log",
    conda:
        "../envs/pacbio.yaml",
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="128000",
        time="24:00:00",
    shell:
        "python {params.scripts}/calc_matrix.py "
        "{input} "
        "{output} "
        "-hamming"

rule fasta_to_levenshtein_distance_matrix:
    input:
        "{base}/aggregated/lineage_clustering/templated/{donor}/{group}_{seq_element}.fasta",
    output:
        "{base}/aggregated/lineage_clustering/templated/{donor}/{group}_{seq_element}.npy",
    log:
        "{base}/logs/cluster_lineages/matrix_calc/templated/{donor}_{group}_{seq_element}.log",
    conda:
        "../envs/pacbio.yaml",
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="128000",
        time="48:00:00",
    shell:
        "python {params.scripts}/calc_matrix.py "
        "{input} "
        "{output} "

def aggregate_cdr3_npy_files(wildcards):
    checkpoint_output = checkpoints.prepare_cdr3_groups_for_distance_evaluation.get(**wildcards).output.groups
    #small_groups = expand("{base}/aggregated/lineage_clustering/cdr3/{donor}/{group}_cdr3.npy",
    #        base=wildcards.base,
    #        donor=wildcards.donor,
    #       group=glob_wildcards(os.path.join(checkpoint_output, "{group}_cdr3.npy")).group)
    large_groups = expand("{base}/aggregated/lineage_clustering/cdr3/{donor}/{group}_cdr3.npy",
            base=wildcards.base,
            donor=wildcards.donor,
           group=glob_wildcards(os.path.join(checkpoint_output, "{group}_cdr3.fasta")).group)
    return large_groups

def aggregate_templated_npy_files(wildcards):
    checkpoint_output = checkpoints.cluster_cdr3s.get(**wildcards).output.clusters
    v_files = expand("{base}/aggregated/lineage_clustering/templated/{donor}/{group}_templated_v.npy",
            base=wildcards.base,
            donor=wildcards.donor,
            group=glob_wildcards(os.path.join(checkpoint_output, "{group}_templated_v.fasta")).group)
    j_files = expand("{base}/aggregated/lineage_clustering/templated/{donor}/{group}_templated_j.npy",
            base=wildcards.base,
            donor=wildcards.donor,
           group=glob_wildcards(os.path.join(checkpoint_output, "{group}_templated_j.fasta")).group)
    return v_files + j_files




# an aggregation over all produced clusters
checkpoint cluster_cdr3s:
    input:
        airr=rules.call_germlines.output.preprocessed,
        npy=aggregate_cdr3_npy_files,
    output:
        airr="{base}/aggregated/lineage_clustering/templated/{donor}/{donor}_unique_vdjs_cdr3_clusters.tsv.gz",
        clusters=directory("{base}/aggregated/lineage_clustering/templated/{donor}/")
    wildcard_constraints:
        donor="TBd[0-9]",
    log:
        "{base}/logs/cluster_cdr3s/{donor}_cluster_cdr3s.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="131000",
        time="24:00:00",
    shell:
        """
        mkdir -p {output.clusters}
        rm {output.clusters}/*

        python {params.scripts}/cluster_sequences_based_on_cdr3_identity_uint8.py \
        {input.airr} \
        -outdir {wildcards.base}/aggregated/lineage_clustering/templated/{wildcards.donor} \
        -matrixdir {wildcards.base}/aggregated/lineage_clustering/cdr3/{wildcards.donor} \
        -samplename {wildcards.donor} \
        2> {log}
        
        """

rule cluster_templated_regions:
    input:
        airr=rules.call_germlines.output.preprocessed,
        cdr3_cluster_table="{base}/aggregated/lineage_clustering/templated/{donor}/{donor}_unique_vdjs_cdr3_clusters.tsv.gz",
        npy=aggregate_templated_npy_files,
    output:
        airr="{base}/aggregated/lineage_clustering/final_lineage_ids/{donor}.tsv.gz",
    log:
        "{base}/logs/cluster_templated/{donor}_cluster_templated.log",
    conda:
        "../envs/pacbio.yaml",
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="131000",
        time="12:00:00",
    shell:
        "python {params.scripts}/cluster_templated_regions_uint8.py "
        "-airr {input.airr} "
        "-cdr3clusters {input.cdr3_cluster_table} "
        "-outdir {wildcards.base}/aggregated/lineage_clustering/final_lineage_ids/ "
        "-matrixdir {wildcards.base}/aggregated/lineage_clustering/templated/{wildcards.donor} "
        "-samplename {wildcards.donor} "
        "2> {log}"

rule create_germline_db:
    input:
        rules.call_germlines.output.germlines,
    output:
        "{base}/aggregated/germline_databases/initial/{donor}_v_germlines.fasta",
    log:
        "{base}/logs/makedb/{donor}_create_germ_db.log",
    conda:
        "../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
    shell:
        "python {params.scripts}/create_germline_database.py "
        "{input} "
        "-outdir {wildcards.base}/aggregated/germline_databases/initial "
        "-samplename {wildcards.donor} "
        "> {log}"


rule blast_to_germline:
    input:
        seqs=rules.cluster_templated_regions.output,
        db=rules.create_germline_db.output,
    output:
        "{base}/aggregated/germline_db_vcall/initial/{donor}.tsv.gz",
    params:
        scripts=config["vdj_scripts"],
    log:
        "{base}/logs/blast_to_germline/{donor}_blast_to_germline.log",
    resources:
        mem_mb="65000",
        time="12:00:00",
    conda:
        "../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/blast_to_germline.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-outdir {wildcards.base}/aggregated/germline_db_vcall/initial "
        "-samplename {wildcards.donor} "
        "2> {log}"


rule polish_germlines:
    input:
        seqs=rules.blast_to_germline.output,
        db=rules.create_germline_db.output,
    output:
        db="{base}/aggregated/germline_databases/final/{donor}_v_germlines_polished.fasta",
    params:
        scripts=config["vdj_scripts"],
        organism=lambda wildcards: samplesheets_vdj[
                samplesheets_vdj.donor == str(wildcards.donor)]["species"].values[0],
        IGDBDIR=config["IGDBDIR"],
    log:
        "{base}/logs/polish_germlines/{donor}_polish_germlines.log",
    conda:
        "../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/polish_germline_db.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-imgt_db {params.IGDBDIR}/fasta/imgt_{params.organism}_ig_v.fasta "
        "-imgt_allele_info {params.IGDBDIR}/internal_data/{params.organism}/{params.organism}.ndm.imgt "
        "-outdir {wildcards.base}/aggregated/germline_databases/final "
        "2> {log}"


rule realign_to_polished_germline:
    input:
        seqs=rules.cluster_templated_regions.output,
        db=rules.polish_germlines.output.db,
    output:
        "{base}/aggregated/germline_db_vcall/final/{donor}.tsv.gz",
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="65000",
        time="12:00:00",
    log:
        "{base}/logs/blast_to_germline/{donor}_blast_to_germline_polished.log",
    conda:
        "../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/blast_to_germline.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-outdir {wildcards.base}/aggregated/germline_db_vcall/final "
        "-samplename {wildcards.donor} "
        "2> {log}"


rule align_v_sequences:
    input:
        seqs=rules.realign_to_polished_germline.output,
        db=rules.polish_germlines.output.db,
    output:
        "{base}/aggregated/vtrees/pseudobulk/{donor}_vmsa.tsv.gz",
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="262000",
        time="2-00:00:00",
    log:
        "{base}/logs/trees/{donor}_pseudobulk_vmsa.log",
    conda:
        "../envs/presto.yaml"
    threads: 20
    shell:
        "python {params.scripts}/align_v_sequences.py "
        "{input.seqs} "
        "-outdir {wildcards.base}/aggregated/vtrees/pseudobulk "
        "-scratchdir {wildcards.base}/aggregated/vtrees/pseudobulk "
        "-samplename {wildcards.donor} "
        "-germline_db {input.db} "
        "-threads {threads} "
        "2> {log}"


rule build_v_trees:
    input:
        rules.align_v_sequences.output,
    output:
        "{base}/aggregated/vtrees/pseudobulk/{donor}_v_trees.tsv",
    params:
        scripts=config["vdj_scripts"],
    log:
        "{base}/logs/trees/{donor}_fasttree.log",
    resources:
        mem_mb="65000",
        time="12:00:00",
    conda:
        "../envs/presto.yaml"
    shell:
        "python {params.scripts}/build_v_trees.py "
        "{input} "
        "-outdir {wildcards.base}/aggregated/vtrees/pseudobulk "
        "-scratchdir {wildcards.base}/aggregated/vtrees/pseudobulk "
        "-samplename {wildcards.donor} "
        "2> {log}"
