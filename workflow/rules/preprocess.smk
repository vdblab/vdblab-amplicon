import os
import traceback

from contextlib import redirect_stderr


include: "common.smk"


envvars:
    "TMPDIR",


min_version("6.0")

#validate(config, os.path.join(workflow.current_basedir, "../config/preprocess.schema.yaml"))

MANIFEST = load_manifest(config["manifest"], None)
SAMPLES = get_samples_from_manifest(MANIFEST)

LOG_PREFIX = "logs/preprocess"


onstart:
    with open("preprocess_config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,


rule all:
    input:
        f"preprocess/{config['pool']}_manifest.tsv",
        f"preprocess/{config['pool']}_missing_or_incomplete.tsv",
        f"preprocess/{config['pool']}_failed_samples.tsv",
        f"preprocess/multiqc_reports/{config['pool']}_dada2_trimming_report_mqc.out",
        f"preprocess/multiqc_reports/{config['pool']}_adapter_contam_report_mqc.out",


rule remove_primers:
    """
    Remove primers from fastq files

    It seems reads can either start with the forward or reverse primer,
    so we provide the -g (--front) arg twice  and --G twice intentionally
    see here:
    https://cutadapt.readthedocs.io/en/stable/guide.html?highlight=reverse%20complement#multiple-adapters

    to convince yourself this is the right move, try the following
    zcat < pool_R1.part_001.fastq.gz | head -n 20 > ./tmp_r_raw.fastq
    zcat < pool_R2.part_001.fastq.gz | head -n 20 > ./tmp_f_raw.fastq

    The following outputs 3 of the 5 reads
    docker run --rm -v $PWD:/data/ ghcr.io/vdblab/cutadapt:3.7 \

    cutadapt \
    -g AYTGGGYDTAAAGNG -G CCGTCAATTYHTTTRAGT \
    -o /data/out_tmp_f -p /data/out_tmp_r \
    --discard-untrimmed \
    /data/tmp_f_raw.fastq /data/tmp_r_raw.fastq

    Where as the following outputs all 5 reads, accounting for those in both orientations
    docker run --rm -v $PWD:/data/ ghcr.io/vdblab/cutadapt:3.7 \
    cutadapt \
    -g AYTGGGYDTAAAGNG -g CCGTCAATTYHTTTRAGT \
    -G AYTGGGYDTAAAGNG -G CCGTCAATTYHTTTRAGT \
    -o /data/out_tmp_f -p /data/out_tmp_r \
    --discard-untrimmed \
    /data/tmp_f_raw.fastq /data/tmp_r_raw.fastq \

    If using this pipeline with samples collected with a different protocol,
    this may need to be modified
    """
    input:
        R1=lambda wildcards: MANIFEST.loc[wildcards.sample, "R1"],
        R2=lambda wildcards: MANIFEST.loc[wildcards.sample, "R2"],
    output:
        R1=temp("preprocess/primers-removed/{sample}_noprimers_R1.fastq.gz"),
        R2=temp("preprocess/primers-removed/{sample}_noprimers_R2.fastq.gz"),
    params:
        primer_F=config["primer_F"],
        primer_R=config["primer_R"],
    container:
        "docker://ghcr.io/vdblab/cutadapt:3.7"
    threads: 4
    resources:
        mem_mb=4 * 1024,
    log:
        o=f"{LOG_PREFIX}/cutadapt_remove_primers_{{sample}}.o",
        e=f"{LOG_PREFIX}/cutadapt_remove_primers_{{sample}}.e",
    shell:
        """
        cutadapt \
            -Z \
            --cores {threads} \
            -g correct={params.primer_F} \
            -g flipped={params.primer_R} \
            -G flipped={params.primer_F} \
            -G correct={params.primer_R} \
            -e 1.5 \
            --minimum-length 150 \
            -o {output.R1} \
            -p {output.R2} \
            --rename='{{header}}  {{adapter_name}}' \
            {input.R1} {input.R2} \
            > {log.o} 2>> {log.e}
        """


rule quality_trim:
    input:
        R1="preprocess/primers-removed/{sample}_noprimers_R1.fastq.gz",
        R2="preprocess/primers-removed/{sample}_noprimers_R2.fastq.gz",
    output:
        R1="preprocess/trimmed/{sample}_trimmed_R1.fastq.gz",
        R2="preprocess/trimmed/{sample}_trimmed_R2.fastq.gz",
        figpath="preprocess/plots/{sample}_quality_profile_mqc.png",
        stats=temp("preprocess/dada2/{sample}_quality_trimming_stats.txt"),
    params:
        trunclen_R1=config["trunclen_R1"],
        trunclen_R2=config["trunclen_R2"],
        figdir_name="dada2/figures/",
        min_asv_len=config["min_asv_len"],
    log:
        e=f"{LOG_PREFIX}/dada2_trim_{{sample}}.e",
        o=f"{LOG_PREFIX}/dada2_trim_{{sample}}.o",
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 2
    resources:
        mem_mb=4 * 1024,
    script:
        "../scripts/preprocess/dada2_trim.R"


rule output_manifest:
    """
    Output a file with sample ids extracted from the oligo files,
    and trimmed R1 and R1 fastq files
    """
    input:
        expand(
            "preprocess/trimmed/{sample}_trimmed_R{dir}.fastq.gz",
            sample=SAMPLES,
            dir=[1, 2],
        ),
    output:
        manifest=f"preprocess/{config['pool']}_manifest.tsv",
        missing=f"preprocess/{config['pool']}_missing_or_incomplete.tsv",
    log:
        e=f"{LOG_PREFIX}/output_manifest.e",
    run:
        with open(log.e, "w") as ef, redirect_stderr(ef):
            try:
                fastq_template = os.path.join(os.getcwd(), "preprocess/trimmed/{sample}_trimmed_R{dir}.fastq.gz")
                write_manifest_and_missing(
                    SAMPLES, fastq_template, output["manifest"], output["missing"]
                )
            except Exception as e:
                traceback.print_exc(file=ef)


rule sample_fastqc_report:
    input:
        R1="preprocess/primers-removed/{sample}_noprimers_R1.fastq.gz",
        R2="preprocess/primers-removed/{sample}_noprimers_R2.fastq.gz",
    output:
        report_R1="preprocess/fastqc_reports/{sample}_noprimers_R1_fastqc.html",
        report_R2="preprocess/fastqc_reports/{sample}_noprimers_R2_fastqc.html",
        zip_R1="preprocess/fastqc_reports/{sample}_noprimers_R1_fastqc.zip",
        zip_R2="preprocess/fastqc_reports/{sample}_noprimers_R2_fastqc.zip",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
    container:
        "docker://staphb/fastqc:0.11.9"
    threads: 2
    resources:
        mem_mb=4 * 1024,
    log:
        e=f"{LOG_PREFIX}/fastqc_{{sample}}.e",
        o=f"{LOG_PREFIX}/fastqc_{{sample}}.o",
    shell:
        """
        fastqc \
            --outdir {params.outdir} \
            --threads {threads} \
            --noextract \
            {input.R1} {input.R2} \
            2> {log.e} > {log.o}
        """


rule get_contam_from_fastqc:
    """
    this prints the overrepresented sequences section with awk and sed,
    then removes all the "No Hit" entries that are likely biologically overrepresented,
    such as the amplicon's primer binding site, non-varying regions between V4-V5, etc
    """
    input:
        zip_R1="preprocess/fastqc_reports/{sample}_noprimers_R1_fastqc.zip",
        zip_R2="preprocess/fastqc_reports/{sample}_noprimers_R2_fastqc.zip",
    output:
        o="preprocess/fastqc_contam/{sample}_fastqc_contaminants.tsv",
    threads: 2
    resources:
        mem_mb=4 * 1024,
    log:
        e=f"{LOG_PREFIX}/get_contam_from_fastqc_{{sample}}.e",
        o=f"{LOG_PREFIX}/get_contam_from_fastqc_{{sample}}.o",
    script:
        "../scripts/preprocess/get_contam_from_fastqc.py"


rule make_dada2_trimming_multiqc_report:
    """
    This is its own rule purely so we can mark the input files as temporary,
    which will make generating unit tests possible
    """
    input:
        trim_reports=expand(
            "preprocess/dada2/{sample}_quality_trimming_stats.txt", sample=SAMPLES
        ),
    output:
        report=f"preprocess/multiqc_reports/{config['pool']}_dada2_trimming_report_mqc.out",
    threads: 1
    resources:
        mem_mb=1 * 1024,
    log:
        e=f"{LOG_PREFIX}/merge_trimming_reports.e",
        o=f"{LOG_PREFIX}/merge_trimming_reports.o",
    shell:
        """
        echo -e "# id: 'trimming'" > {output.report}
        echo -e "# plot_type: 'table'" >> {output.report}
        echo -e "# section_name: 'DADA2 Quality Trimming'" >> {output.report}
        echo -e "sample\tprefilter\tpostfilter\tpct_loss" >> {output.report}
        cat {input.trim_reports} | grep -v "reads.in" >> {output.report} 2>> {log.e}
        """


rule make_fastqc_contam_multiqc_report:
    """
    This is its own rule purely so we can mark the input files as temporary,
    which will make generating unit tests possible
    """
    input:
        fastqc_reports=expand(
            "preprocess/fastqc_contam/{sample}_fastqc_contaminants.tsv", sample=SAMPLES
        ),
    output:
        report=f"preprocess/multiqc_reports/{config['pool']}_adapter_contam_report_mqc.out",
    threads: 1
    resources:
        mem_mb=1 * 1024,
    log:
        e=f"{LOG_PREFIX}/make_fastqc_contam_multiqc_report.e",
        o=f"{LOG_PREFIX}/make_fastqc_contam_multiqc_report.o",
    shell:
        """
        echo -e "# id: 'fastqc_adapters'" > {output.report}
        echo -e "# plot_type: 'table'" >> {output.report}
        echo -e "# section_name: 'FastQC-detected Adapter Contamination'" \
            >> {output.report}
        head -n 1 \
            {input.fastqc_reports[0]} >> {output.report} 2>> {log.e}
        tail -n +2 -q \
            {input.fastqc_reports} >> {output.report} 2>> {log.e}
        """


rule detect_failed_samples:
    input:
        trim_report=f"preprocess/multiqc_reports/{config['pool']}_dada2_trimming_report_mqc.out",
        adapter_contam_report=f"preprocess/multiqc_reports/{config['pool']}_adapter_contam_report_mqc.out",
    output:
        out=f"preprocess/{config['pool']}_failed_samples.tsv",
    params:
        min_retained=config["min_retained"],
        min_read_pairs=config["min_read_pairs"],
        max_perc_adapter=config["max_perc_adapter"],
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 1
    resources:
        mem_mb=1 * 1024,
    log:
        e=f"{LOG_PREFIX}/detect_failed_samples.e",
    script:
        "../scripts/preprocess/detect_failed_samples.R"
