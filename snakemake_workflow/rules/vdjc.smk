samplesheet = samplesheets
import os

rule igblast:
    input:
        "{base}/immcantation/{sample_uid}/filtered_contig.fasta"
    output:
        "{base}/igblast/{sample_uid}_igblast.tsv",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        organism=lambda wildcards: samplesheet_lookup(wildcards.sample_uid, "species"),
        IGDBDIR=config["IGDBDIR"],
    shell:
        """
        
        ml load system libuv
        wdir=$(dirname {input[0]})
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
        -out {output}
        """


def aggregate_input(wildcards):
    samplelist = [
        "{}/igblast/{}_igblast.tsv".format(config["outs_basedir"], bc)
        for bc in samplesheets[samplesheets.donor == wildcards.sample].index
    ]
    return samplelist

rule combine_samples:
    input:
        aggregate_input,
    output:
        "{base}/combined/{sample}_combined.tsv.gz",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    log:
        "{base}/logs/{sample}_combined.log",
    params:
        scripts=config["scripts"],
        samplesheet=config["samplesheets"],
    shell:
        """
         python {params.scripts}/combine_igblast_output.py \
        {input} \
        -samplesheet {params.samplesheet} \
         -outdir {wildcards.base}/combined \
         -name {wildcards.sample}
        """


rule filter_vdj:
    input:
        "{base}/combined/{sample}_combined.tsv.gz",
    output:
        "{base}/vdj/{sample}_combined_vdj_filtered.tsv.gz",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    log:
        "{base}/logs/{sample}_filter_vdj.log",
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/filter_igblast_output.py "
        "{input} "
        "-outdir {wildcards.base}/vdj "
        "--verbose "
        "> {log} "
	"2> {log}"

rule filter_vdj_passthru:
    input:
        "{base}/combined/{sample}_combined.tsv.gz",
    output:
        "{base}/vdj/{sample}_combined_vdj.tsv.gz",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    log:
        "{base}/logs/{sample}_filter_vdj.log",
    params:
        scripts=config["scripts"],
    shell:
        "cp {input} {output}"


rule annotate_constant_region:
    input:
        "{base}/vdj/{sample}_combined_vdj.tsv.gz",
    output:
        "{base}/vdjc/{sample}_combined_vdjc.tsv.gz",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        ighc_db=lambda wildcards: config["ighc_db"][
            samplesheets[samplesheets.donor == str(wildcards.sample)]["species"].values[
                0
            ]
        ],
        scripts=config["scripts"],
    log:
        "{base}/logs/{sample}_annotate_constant_region.log",
    shell:
        "python {params.scripts}/blast_constant_region.py "
        "{input} "
        "-ighc_db {params.ighc_db} "
        "-outdir {wildcards.base}/vdjc "
        "> {log}"

rule get_creads:
    input:
        "{base}/vdjc/{sample}_combined_vdjc.tsv.gz",
    output:
        "{base}/vdjc/{sample}_combined_creads.fasta"
    log:
        "{base}/logs/{sample}_combined_creads.log",
    conda:
       os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/c_reads_2_fasta.py "
        "{input} "
        "{output} "
        "2> {log}"

rule call_germlines:
    input:
        "{base}/vdjc/{sample}_combined_vdjc.tsv.gz",
    output:
        germlines="{base}/grmlin/{sample}_combined_vdjc_v_germlines.tsv",
	preprocessed="{base}/grmlin/{sample}_combined_vdjc_preprocessed.tsv.gz",
    log:
        "{base}/logs/{sample}_grmlin.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        organism=lambda wildcards: samplesheets[
            samplesheets.donor == str(wildcards.sample)
        ]["species"].values[0],
        grmlin=config["path_to_grmlin"],
        IGDBDIR=config["IGDBDIR"],
    shell:
        "{params.grmlin}/grmlin "
        "{input} "
        "-annotate {params.IGDBDIR}/database/imgt_{params.organism}_ig_v "
        "-outdir {wildcards.base}/grmlin "
        "--verbose "
        "> {log}"

rule cluster_lineages:
    input:
        "{base}/grmlin/{sample}_combined_vdjc_preprocessed.tsv.gz",
    output:
        "{base}/lineages/{sample}_combined_vdjc_lineage_ids.tsv.gz",
    log:
        "{base}/logs/{sample}_cluster_lineages.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/cluster_lineages_small_dataset.py "
        "{input} "
        "-outdir {wildcards.base}/lineages "
        "2> {log}"


rule create_germline_db:
    input:
        "{base}/grmlin/{sample}_combined_vdjc_v_germlines.tsv",
    output:
        "{base}/germline_databases/{sample}_v_germlines.fasta",
    log:
        "{base}/logs/{sample}_create_germ_db.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/create_germline_database.py "
        "{input} "
        "-outdir {wildcards.base}/germline_databases "
        "-samplename {wildcards.sample} "
        "> {log}"


rule blast_to_germline:
    input:
        seqs="{base}/lineages/{sample}_combined_vdjc_lineage_ids.tsv.gz",
        db="{base}/germline_databases/{sample}_v_germlines.fasta",
    output:
        "{base}/vseq_aligned_to_germline/{sample}_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz",
    params:
        scripts=config["scripts"],
    log:
        "{base}/logs/{sample}_blast_to_germline.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    shell:
        "python {params.scripts}/blast_to_germline.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-outdir {wildcards.base}/vseq_aligned_to_germline "
        "2> {log}"


rule polish_germlines:
    input:
        seqs="{base}/vseq_aligned_to_germline/{sample}_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz",
        db="{base}/germline_databases/{sample}_v_germlines.fasta",
    output:
        db="{base}/polished_germline_databases/{sample}_v_germlines_polished.fasta",
    params:
        scripts=config["scripts"],
        organism=lambda wildcards: samplesheets[
            samplesheets.donor == str(wildcards.sample)
        ]["species"].values[0],
        IGDBDIR=config["IGDBDIR"],
    log:
        "{base}/logs/{sample}_polish_germlines.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    shell:
        "python {params.scripts}/polish_germline_db.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-imgt_db {params.IGDBDIR}/fasta/imgt_{params.organism}_ig_v.fasta "
        "-imgt_allele_info {params.IGDBDIR}/internal_data/{params.organism}/{params.organism}.ndm.imgt "
        "-outdir {wildcards.base}/polished_germline_databases "
        "2> {log}"


rule realign_to_polished_germline:
    input:
        seqs="{base}/lineages/{sample}_combined_vdjc_lineage_ids.tsv.gz",
        db="{base}/polished_germline_databases/{sample}_v_germlines_polished.fasta",
    output:
        "{base}/vseq_aligned_to_germline_polished/{sample}_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz",
    params:
        scripts=config["scripts"],
    log:
        "{base}/logs/{sample}_blast_to_germline_polished.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    shell:
        "python {params.scripts}/blast_to_germline.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-outdir {wildcards.base}/vseq_aligned_to_germline_polished "
        "2> {log}"


rule align_v_sequences:
    input:
        seqs="{base}/vseq_aligned_to_germline_polished/{sample}_combined_vdjc_lineage_ids_vseq_germ_blast.tsv.gz",
        db="{base}/polished_germline_databases/{sample}_v_germlines_polished.fasta",
    output:
        "{base}/trees/{sample}_vmsa.tsv.gz",
    params:
        scripts=config["scripts"],
    log:
        "{base}/logs/{sample}_vmsa.log",
    conda:
        os.path.join(workflow.basedir, "envs/presto.yaml")
    shell:
        "python {params.scripts}/align_v_sequences.py "
        "{input.seqs} "
        "-outdir {wildcards.base}/trees "
        "-scratchdir {wildcards.base}/trees "
        "-samplename {wildcards.sample} "
        "-germline_db {input.db} "
        "2> {log}"


rule build_v_trees:
    input:
        "{base}/trees/{sample}_vmsa.tsv.gz",
    output:
        "{base}/trees/{sample}_v_trees.tsv",
    params:
        scripts=config["scripts"],
    log:
        "{base}/logs/{sample}_fasttree.log",
    conda:
        os.path.join(workflow.basedir, "envs/pacbio.yaml")
    shell:
        "python {params.scripts}/build_v_trees.py "
        "{input} "
        "-outdir {wildcards.base}/trees "
        "-scratchdir {wildcards.base}/trees "
        "-samplename {wildcards.sample} "
        "2> {log}"
