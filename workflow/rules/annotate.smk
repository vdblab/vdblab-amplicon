include: "common.smk"


envvars:
    "TMPDIR",


min_version("6.0")

validate(config, os.path.join(workflow.basedir, "../config/annotate.schema.yaml"))

config["pipeline_version"] = get_pipeline_version()

LOG_PREFIX = "logs/annotate"


onstart:
    with open("annotate_config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,


rule all:
    input:
        f"annotate/{config['pool']}_asv_annotated.txt",
        f"annotate/{config['pool']}_asv_blast_out.tsv",
        f"annotate/multiqc_reports/{config['pool']}_host_contamination_mqc.out",


rule detect_contamination:
    """
    Parameters from LOTUS2, even though
    https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01012-1
    says default params.
    """
    input:
        asv_fasta=config["asv_fasta"],
        ref=lambda wildcards: config[f"{wildcards.ref}_ref"],
    output:
        hits=f"annotate/{config['pool']}_{{ref}}_hits",
    container:
        "docker://evolbioinfo/minimap2:v2.24"
    resources:
        mem_mb=lambda wildcards: (
            18 * 1024 if wildcards.ref in ["human", "mouse"] else 2 * 1024
        ),
    log:
        e=f"{LOG_PREFIX}/search_host_content_{{ref}}.e",
    threads: 16
    shell:
        """
        minimap2 \
            -x sr \
            --sr \
            -u both \
            --secondary=no \
            -N 30 \
            -c \
            -t {threads} \
            -o {output.hits} \
            {input.ref} {input.asv_fasta} \
            2>> {log.e}
        """


rule make_contamination_multiqc_table:
    """
    make a multiqc output table
    In the future maybe it should output the contamination table as well;
    the multiqc header could be added to empty files and the minimap results
    appended, rather than using the -o syntax for output
    """
    input:
        **(
            {"adapter_hits": f"annotate/{config['pool']}_adapter_hits"}
            if config["skip_contam_search"]
            else {
                "adapter_hits": f"annotate/{config['pool']}_mouse_hits",
                "human_hits": f"annotate/{config['pool']}_human_hits",
                "mouse_hits": f"annotate/{config['pool']}_adapter_hits",
            }
        ),
    output:
        report=f"annotate/multiqc_reports/{config['pool']}_host_contamination_mqc.out",
    resources:
        mem_mb=1 * 1024,
    threads: 1
    log:
        e=f"{LOG_PREFIX}/make_contamination_multiqc_table.e",
    script:
        "../scripts/annotate/make_contamination_multiqc_table.py"


rule legacy_blast_annotate:
    input:
        asv_fasta=config["asv_fasta"],
    output:
        blast=f"annotate/{config['pool']}_asv_blast_out.tsv",
    resources:
        mem_mb=4 * 1024,
    threads: 8
    container:
        "docker://ghcr.io/vdblab/micro16s-blast:2.13.0a"
    log:
        e=f"{LOG_PREFIX}/legacy_blast_annotate.e",
        o=f"{LOG_PREFIX}/legacy_blast_annotate.o",
    shell:
        """
        blastn \
            -num_threads {threads} \
            -db /db/16S_ribosomal_RNA \
            -max_target_seqs 500 \
            -query {input.asv_fasta} \
            -outfmt \"6 qseqid staxids saccver stitle qlen length nident pident bitscore evalue score\" \
            -out {output.blast} \
            >> {log.o} 2> {log.e}
        """


rule parse_blast_annotations:
    input:
        blast=f"annotate/{config['pool']}_asv_blast_out.tsv",
        annotation=config["anno_db"],
    output:
        annotations=f"annotate/{config['pool']}_asv_annotated.txt",
        passed=f"annotate/{config['pool']}_blast_passed.txt",
        not_passed=f"annotate/{config['pool']}_blast_not_passed.txt",
        detailed=f"annotate/{config['pool']}_blast_passed_not_passed.txt",
    threads: 1
    resources:
        mem_mb=2 * 1024,
    container:
        "docker://rocker/tidyverse:4.1"
    log:
        e=f"{LOG_PREFIX}/parse_blast_hits.e",
        o=f"{LOG_PREFIX}/parse_blast_hits.o",
    script:
        "../scripts/annotate/parse_blast_hits.R"
