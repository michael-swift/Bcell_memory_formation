smplesheet = samplesheets


include: "vdjc.smk"


rule call_cells:
    input:
        VDJ_file=rules.realign_to_polished_germline.output,
   output:
        called_cells="{base}/aggregated/cell_calls/{donor}_called_cells.tsv.gz",
        ambient_rna="{base}/aggregated/cell_calls/{donor}_ambient_vdjs.tsv.gz",
    log:
        "{base}/logs/{donor}_call_cells.log",
    conda:
        "../envs/pacbio.yaml"
    resources:
        mem_mb="65000",
    params:
        scripts=config["vdj_scripts"],
        whitelist=config["whitelist"],
    shell:
        "python {params.scripts}/call_cells.py {input.VDJ_file} "
        "-whitelist {params.whitelist} "
        "-outdir {wildcards.base}/aggregated/cell_calls "
        "-outname {wildcards.donor} "
        "2> {log}"


rule combine_cells_with_vdj_annotations:
    input:
        vdj_info=rules.realign_to_polished_germline.output,
        cell_calls=rules.call_cells.output.called_cells,
    output:
        "{base}/aggregated/cell_calls/{donor}_called_cells_vdj_annotated.tsv.gz",
    log:
        "{base}/logs/{donor}_combined_cell_call_and_annotation.log",
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="262000",
    shell:
        "python {params.scripts}/combine_called_cell_with_vseq_data.py "
        "-cell_calls {input.cell_calls} "
        "-vdj_info {input.vdj_info} "
        "-outdir {wildcards.base}/aggregated/cell_calls "
        "-outname {wildcards.donor}_called_cells_vdj_annotated "
        "2> {log}"

rule add_sample_info:
    input:
        rules.combine_cells_with_vdj_annotations.output,
    output:
        "{base}/aggregated/cell_calls/{donor}_called_cells_vdj_annotated_extended.tsv.gz",
    log:
        "{base}/logs/{donor}_add_sample_info.log",
    params:
        scripts=config["vdj_scripts"],
        samplesheet=config["samplesheets"][0],
    resources:
        mem_mb="64000",
    shell:
        "python {params.scripts}/add_sample_info.py "
        "{input} "
        "-samplesheet {params.samplesheet} "
        "-output {output} "
        "> {log}"

rule align_cell_v_sequences:
    input:
        seqs=rules.combine_cells_with_vdj_annotations.output,
        db=rules.polish_germlines.output.db,
    output:
        msa="{base}/aggregated/vtrees/cells/{donor}_vmsa.tsv.gz",
        scratch=directory("{base}/aggregated/vtrees/cells/{donor}_scratch"),
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="64000",
        time="1-00:00:00",
        threads=20,
    log:
        "{base}/logs/trees/{donor}_cell_vmsa.log",
    conda:
        "../envs/presto.yaml"
    threads: 20
    shell:
        """
        mkdir -p {wildcards.base}/aggregated/vtrees/cells/{wildcards.donor}_scratch

        python {params.scripts}/align_v_sequences.py \
        {input.seqs} \
        -outdir {wildcards.base}/aggregated/vtrees/cells \
        -scratchdir {wildcards.base}/aggregated/vtrees/cells/{wildcards.donor}_scratch \
        -samplename {wildcards.donor} \
        -germline_db {input.db} \
        -threads {threads} \
        2> {log}
        """


rule build_cell_v_trees:
    input:
        rules.align_cell_v_sequences.output.msa,
    output:
        "{base}/aggregated/vtrees/cells/{donor}_v_trees.tsv",
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
        "-outdir {wildcards.base}/aggregated/vtrees/cells "
        "-scratchdir {wildcards.base}/aggregated/vtrees/cells/{wildcards.donor}_scratch "
        "-samplename {wildcards.donor} "
        "2> {log}"

def fetch_all_donor_cell_calls(wildcards):
    files = [
        "{}/aggregated/cell_calls/"
        "{}_called_cells_vdj_annotated_extended.tsv.gz".format(
            wildcards.base, donor
        )
        for donor in samplesheets_vdj.donor.unique()
        ]
    return files

rule annotate_likely_cross_contaminants:
    input:
        fetch_all_donor_cell_calls,
    output:
        tsv="{base}/all_vdj_cell_calls_IGH.tsv.gz",
        figures=directory("{base}/figures/cross_contamination_stats"),
    params:
        scripts=config["vdj_scripts"],
    log:
        "{base}/logs/remove_cross_contaminating_vdjs_IGH.log"
    resources:
        mem_meb="65000",
    conda:
        "../envs/scanpy.py"
    shell:
        "python {params.scripts}/annotate_cross-contaminating_barcodes.py "
        "-input_paths {input} "
        "-outname all_vdj_cell_calls "
        "-outdir {wildcards.base} "
        "-figure_outdir {output.figures} "
        "locus IGH"
        "2> {log}"
