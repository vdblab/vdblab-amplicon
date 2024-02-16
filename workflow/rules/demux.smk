import os
import pandas as pd
import traceback

from contextlib import redirect_stderr
from snakemake.utils import validate


include: "common.smk"


envvars:
    "TMPDIR",


wildcard_constraints:
    dir="1|2",


min_version("6.0")


# validate(config, os.path.join(workflow.current_basedir, "../../config/demux.schema.yaml"))


BARCODES = extract_barcodes(config["oligos"])
SAMPLES = BARCODES["sample_id"].tolist()
SHARDS = make_shard_names(config["nshards"])
NLIBS = [i for i, x in enumerate(config["R1"])]
LOG_PREFIX = "logs/demux"


onstart:
    with open("demux_config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,


rule all:
    input:
        f"demux/{config['pool']}_manifest.tsv",
        f"demux/{config['pool']}_missing_or_incomplete.tsv",
        expand(
            f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R1_fastqc.html",
            lib=NLIBS,
        ),
        expand(
            f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R2_fastqc.html",
            lib=NLIBS,
        ),
        expand("demux/fastq/{sample}_R{dir}.fastq.gz", sample=SAMPLES, dir=[1, 2]),


rule concat_fastqs:
    """
    Combine fastq files, and ensure consistent file naming
    """
    input:
        R1=config["R1"],
        R2=config["R2"],
    output:
        R1=temp(f"{config['pool']}_R1.fastq.gz"),
        R2=temp(f"{config['pool']}_R2.fastq.gz"),
    params:
        n_files=len(config["R1"]),
    log:
        e=f"{LOG_PREFIX}/fetch_and_cat_fastqs.e",
        o=f"{LOG_PREFIX}/fetch_and_cat_fastqs.o",
    threads: 1
    resources:
        mem_mb=2 * 1024,
    shell:
        """
        cat {input.R1} > {output.R1} 2> {log.e}
        cat {input.R2} > {output.R2} 2>> {log.e}
        """


rule generate_pool_fastqc_report:
    """ we generate per fastq fastqc reports.
    the symlinking is to ensure consitent names
    """
    input:
        R1=lambda wildcards: config["R1"][int(wildcards.lib)],
        R2=lambda wildcards: config["R2"][int(wildcards.lib)],
    output:
        rep_R1=f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R1_fastqc.html",
        rep_R2=f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R2_fastqc.html",
        zip_R1=f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R1_fastqc.zip",
        zip_R2=f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R2_fastqc.zip",
        R1=temp(f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R1.fq.gz"),
        R2=temp(f"demux/fastqc_reports/{config['pool']}_lib{{lib}}_R2.fq.gz"),
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        base=lambda wildcards, output: os.path.basename(output[0]).replace(
            "1_fastqc.html", ""
        ),
    container:
        "docker://staphb/fastqc:0.11.9"
    # fastqc uses one thread per input file, so increasing this won't help
    threads: 2
    resources:
        mem_mb=4 * 1024,
    log:
        o=f"{LOG_PREFIX}/pool_fastqc_lib{{lib}}.o",
        e=f"{LOG_PREFIX}/pool_fastqc_lib{{lib}}.e",
    shell:
        """
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        fastqc \
            --outdir {params.outdir} \
            --threads {threads} \
            --noextract \
            {output.R1} {output.R2} \
            2> {log.e} > {log.o}
        ls
        """


module utils:
    config:
        config
    snakefile:
        "utils.smk"


use rule split_fastq from utils as split_fastq with:
    input:
        R1=f"{config['pool']}_R1.fastq.gz",
        R2=f"{config['pool']}_R2.fastq.gz",
    output:
        R1=temp(
            expand(
                f"split_fastq/{config['pool']}_R1.part_{{shard}}.fastq.gz",
                shard=SHARDS,
            )
        ),
        R2=temp(
            expand(
                f"split_fastq/{config['pool']}_R2.part_{{shard}}.fastq.gz",
                shard=SHARDS,
            )
        ),
    params:
        nshards=config["nshards"],
        outdir=lambda w, output: directory(os.path.dirname(output.R1[0])),
    log:
        o=f"{LOG_PREFIX}/split_fastq.o",
        e=f"{LOG_PREFIX}/split_fastq.e",


rule make_barcodes_fasta:
    output:
        fasta=f"demux/{config['pool']}_barcodes_R{{dir}}.fasta",
    log:
        e=f"{LOG_PREFIX}/make_barcodes_fasta_{{dir}}.e",
    run:
        with open(log.e, "w") as ef, redirect_stderr(ef):
            try:
                with open(output["fasta"], "w") as out_f:
                    for sample_id, f, r in BARCODES.itertuples(index=False):
                        out_f.write(f">{sample_id}\n")

                        if wildcards.dir == "1":
                            out_f.write(f"{f}\n")
                        else:
                            out_f.write(f"{r}\n")

            except Exception as e:
                traceback.print_exc(file=ef)


rule demux_round_one:
    """
    Demultiplexing round one: check forward reads

    You may be tempted to try to demux and deprimer in the same step;
    Nick tried that but ended with many fewer reads than the two-step process.
    but he may have missed something...

    see https://github.com/benjjneb/dada2/issues/159: need to set a minimum length.
    """
    input:
        barcodes=f"demux/{config['pool']}_barcodes_R1.fasta",
        R1=(
            f"split_fastq/{config['pool']}_R1.part_{{shard}}.fastq.gz"
            if config["nshards"] > 1
            else f"{config['pool']}_R1.fastq.gz"
        ),
        R2=(
            f"split_fastq/{config['pool']}_R2.part_{{shard}}.fastq.gz"
            if config["nshards"] > 1
            else f"{config['pool']}_R2.fastq.gz"
        ),
    output:
        fastq=temp(
            expand(
                "demux/raw/{sample}_R{dir}_round1_s{{shard}}.fastq.gz",
                sample=SAMPLES,
                dir=[1, 2],
            )
        ),
        orphans_R1=temp("demux/raw/orphans_R1_round1_s{shard}.fastq.gz"),
        orphans_R2=temp("demux/raw/orphans_R2_round1_s{shard}.fastq.gz"),
    params:
        # curly-enclosed "name" is imputed by cutadapt itself
        template_R1=lambda wildcards: f"demux/raw/{{name}}_R1_round1_s{wildcards.shard}.fastq.gz",
        template_R2=lambda wildcards: f"demux/raw/{{name}}_R2_round1_s{wildcards.shard}.fastq.gz",
    log:
        o=f"{LOG_PREFIX}/cutadapt_demux_round1_s{{shard}}.o",
        e=f"{LOG_PREFIX}/cutadapt_demux_round1_s{{shard}}.e",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=6 * 60,
    container:
        "docker://ghcr.io/vdblab/cutadapt:3.7"
    shell:
        """
        cutadapt \
            -Z \
            --cores {threads} \
            -e 2 \
            --no-indels \
            --minimum-length 150 \
            -g file:{input.barcodes} \
            -o {params.template_R1} \
            -p {params.template_R2} \
            --untrimmed-output {output.orphans_R1} \
            --untrimmed-paired-output {output.orphans_R2} \
            {input.R1} {input.R2} \
            > {log.o} 2> {log.e}
        """


rule demux_round_two:
    """
    Demultiplexing round two: check reverse reads barcodes among the first round failures

    We provide the inverse orientation of the R1 and R2 input files
    and name the output the same way, so that correct orientation names are preserved

    This doesn't actually separate reads from biological orientation,
    it just gets around issue with detecting barcodes with issues in the forward direction

    This is important because the R1 and R2 reads have very different error profiles,
    so we need to keep the sequencing orientation rather than the biological orientation
    """
    input:
        barcodes=f"demux/{config['pool']}_barcodes_R2.fasta",
        R1="demux/raw/orphans_R1_round1_s{shard}.fastq.gz",
        R2="demux/raw/orphans_R2_round1_s{shard}.fastq.gz",
    output:
        fastq=temp(
            expand(
                "demux/raw/{sample}_R{dir}_round2_s{{shard}}.fastq.gz",
                sample=SAMPLES,
                dir=[1, 2],
            )
        ),
        orphans_R1=temp("demux/raw/orphans_R1_round2_s{shard}.fastq.gz"),
        orphans_R2=temp("demux/raw/orphans_R2_round2_s{shard}.fastq.gz"),
    params:
        # curly-enclosed "name" is imputed by cutadapt itself
        template_R1=lambda wildcards: f"demux/raw/{{name}}_R1_round2_s{wildcards.shard}.fastq.gz",
        template_R2=lambda wildcards: f"demux/raw/{{name}}_R2_round2_s{wildcards.shard}.fastq.gz",
    log:
        o=f"{LOG_PREFIX}/cutadapt_demux_round2_s{{shard}}.o",
        e=f"{LOG_PREFIX}/cutadapt_demux_round2_s{{shard}}.e",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=6 * 60,
    container:
        "docker://ghcr.io/vdblab/cutadapt:3.7"
    shell:
        """
        cutadapt \
            -Z \
            --cores {threads} \
            -e 2 \
            --no-indels \
            -g file:{input.barcodes} \
            -o {params.template_R2} \
            -p {params.template_R1} \
            --minimum-length 150 \
            --untrimmed-output {output.orphans_R2} \
            --untrimmed-paired-output {output.orphans_R1} \
            {input.R2} {input.R1} \
            > {log.o} 2> {log.e}
        """


rule merge_shards:
    """
    Merge fastqs for shards and demultiplexing rounds
    """
    input:
        fastqs=lambda wildcards: (
            expand(
                "demux/raw/{{sample}}_R{{dir}}_round{round}_s{shard}.fastq.gz",
                    round=[1, 2],
                    shard=SHARDS,
                )
            if wildcards.sample != "orphans"
            else expand(
                "demux/raw/{{sample}}_R{{dir}}_round2_s{shard}.fastq.gz",
                shard=SHARDS,
            )
        ),
    output:
        fastqs="demux/fastq/{sample}_R{dir}.fastq.gz",
    threads: 1
    resources:
        mem_mb=1 * 1024,
    log:
        e=f"{LOG_PREFIX}/merge_shards_{{sample}}_R{{dir}}.e",
    shell:
        """
        cat {input} > {output} 2> {log.e}
        """


rule output_manifest:
    """
    Output a file with sample ids extracted from the oligo files,
    and demultiplexed R1 and R1 fastq files
    """
    input:
        expand(
            "demux/fastq/{sample}_R{dir}.fastq.gz",
            sample=SAMPLES + ["orphans"],
            dir=[1, 2],
        ),
    output:
        manifest=f"demux/{config['pool']}_manifest.tsv",
        missing=f"demux/{config['pool']}_missing_or_incomplete.tsv",
    log:
        e=f"{LOG_PREFIX}/output_manifest.e",
    run:
        with open(log.e, "w") as ef, redirect_stderr(ef):
            try:
                sample_ids = SAMPLES + ["orphans"]
                fastq_template = os.path.join(
                    os.getcwd(), "demux/fastq/{sample}_R{dir}.fastq.gz"
                )
                write_manifest_and_missing(
                    sample_ids, fastq_template, output["manifest"], output["missing"]
                )
            except Exception as e:
                traceback.print_exc(file=ef)
