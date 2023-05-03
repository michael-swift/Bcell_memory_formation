import os


rule get_genome_fasta_and_gtf:
    """ Get genome and GTF from ensembl """
    output:
        fasta="%s/resources/star/%s.ERCC.fa"
        % (workflow.basedir, config["fasta_name"][config["species"]]),
        gtf="%s/resources/star/%s.ERCC.gtf"
        % (workflow.basedir, config["fasta_name"][config["species"]]),
    threads: 1
    params:
        name="wget_files",
        partition=config["partition"],
        fasta_url=config["fasta_url"][config["species"]],
        fasta_name=config["fasta_name"][config["species"]],
        gtf_url=config["gtf_url"][config["species"]],
        gtf_name=config["gtf_name"][config["species"]],
        include_ERCCs="true" if config["include_ERCCs"] else "false",
    resources:
        mem_mb=5000,
    shell:
        "fa=resources/star/{params.fasta_name}.fa.gz && "
        "gtf=resources/star/{params.gtf_name}.gtf.gz && "
        "wget {params.fasta_url} -O $fa && "
        "wget {params.gtf_url} -O $gtf && "
        "zcat $fa > {output.fasta} && "
        "zcat $gtf > {output.gtf} && "
        "if [ {params.include_ERCCs} = true ]; then "
        "cat resources/star/ERCC92.fa >> {output.fasta} && "
        "cat resources/star/ERCC92.gtf >> {output.gtf}; fi"


rule star_genome_generate:
    """ Build the STAR genome
        Notes:
        * --genomeSAsparseD is used to reduce genome memory consumption
            during mapping
    """
    input:
        fasta=rules.get_genome_fasta_and_gtf.output.fasta,
        gtf=rules.get_genome_fasta_and_gtf.output.gtf,
        sjadditional="{}/resources/star/IGHC_IGHJ_splices_{}.txt".format(
            workflow.basedir, config["species"]
        ),
    output:
        "%s/resources/star/star_genome_%s/Genome"
        % (workflow.basedir, config["species"]),
    threads: 12
    params:
        name="star_genome_gen",
        star_read_len=int(config["read_length"]) - 1,
    resources:
        mem_mb=60000,
        partition="quake,owners",
    conda:
        config["workflow_dir"] + "/envs/single-cell.yaml"
    shell:
        "genomedir=$(dirname {output}) && "
        "mkdir -p $genomedir && "
        "cd $(dirname $genomedir) && "
        "STAR --runMode genomeGenerate "
        "--genomeDir $genomedir "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "--sjdbOverhang {params.star_read_len} "
        "--sjdbFileChrStartEnd {input.sjadditional} "
        "--runThreadN {threads}"


rule pull_TICA_data:
    """ Get TICA data from tissueimmunecellatlas.org """
    output:
        output_dir=directory("{base}/downloads/"),
        global_counts_data="{base}/downloads/CountAdded_PIP_global_object_for_cellxgene.h5ad",
        GEX_BCR_data="{base}/downloads/TICA_B_BCR.h5ad",
    threads: 1
    params:
        name="wget_files",
        partition=config["partition"],
        global_url=config["TICA_urls"]["global"],
        bcell_url=config["TICA_urls"]["b_cell"],
    resources:
        mem_mb=5000,
    shell:
        "mkdir -p {output[0]} && "
        "cd {output[0]} && "
        "wget {params.global_url} && "
        "wget {params.bcell_url}"


rule add_TICA_metadata:
    """ add tissue column for integration with our data"""
    input:
        rules.pull_TICA_data.output.global_counts_data,
    output:
        mod_global_counts_data="{base}/downloads/mod_PIP_global_object_for_cellxgene.h5ad.gz",
    threads: 1
    params:
        name="modify_TICA",
        partition=config["partition"],
    resources:
        mem_mb=32000,
    run:
        import scanpy as sc

        TICA_to_Bursa = TICA_to_Bursa = {
            "BLD": "PB",
            "BMA": "BM",
            "SPL": "SP",
            "ILN": "LN",
            "MLN": "LN",
        }
        adata = sc.read_h5ad(input[0])
        adata.obs.loc[:, "tissue"] = adata.obs.Organ.map(
            lambda x: TICA_to_Bursa.get(x, x)
        )
        adata.obs.loc[:, "sample_uid"] = (
            adata.obs.Organ.astype(str) + "_" + adata.obs.Donor.astype(str)
        )
        adata.obs.loc[:, "donor"] = adata.obs.Donor
        # these are not real cellbender counts rn
        adata.X = adata.layers["counts"]
        adata.layers["cellbender_counts"] = adata.X
        adata.write_h5ad(output[0], compression="gzip")


rule star_solo_vdj:
    """ Map reads to genome using STAR
        Notes:
        * Expects gzipped fastqs
        * --twopassMode for better splice alignments
        * --outSAMmapqUnique 60 for any downstream variant analysis
        * --outFilterMismatchNmax disabled by high value
        * --outFilterMismatchNoverReadLmax set higher to account for
            somatic hypermutation
    """
    input:
        rules.star_genome_generate.output,
        get_r1_fqgz_using_wildcards,
        get_r2_fqgz_using_wildcards,
    output:
        "{base}/per_sample/star_solo_vdj/{sample_uid_vdj}/Aligned.out.bam",
    params:
        name="star",
        partition=config["partition"],
        cr_whitelist=config["cr_whitelist"],
    resources:
        mem_mb=60000,
    threads: 12
    conda:
        config["workflow_dir"] + "/envs/single-cell.yaml"
    shell:
        "wdir=$(dirname {output}) && "
        "echo $wdir && "
        "mkdir -p $wdir && "
        "STAR "
        "--genomeDir $(dirname {input[0]}) "
        "--readFilesIn {input[1]} {input[2]} "
        "--readFilesCommand gunzip -c "
        "--outSAMmapqUnique 60 "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverReadLmax 0.1 "
        "--twopassMode Basic "
        "--runThreadN {threads} "
        "--outSAMtype BAM Unsorted "
        "--outFileNamePrefix $wdir/"
        " --soloBarcodeMate 1"
        " --clip5pNbases 39 0"
        " --soloType CB_UMI_Simple"
        " --soloCBstart 1"
        " --soloCBlen 16"
        " --soloUMIstart 17"
        " --soloUMIlen 10"
        " --soloFeatures SJ"
        " --soloCBwhitelist {params.cr_whitelist}"
