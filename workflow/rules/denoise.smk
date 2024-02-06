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

validate(config, os.path.join(workflow.basedir, "../../config/denoise.schema.yaml"))

assert config["pooling"] == "none", f"unpooled workflow selected but config specifies pooling={config['pooling']}"

MANIFEST = load_manifest(config["manifest"], None)
SAMPLES = get_samples_from_manifest(MANIFEST)

LOG_PREFIX = "logs/denoise"


onstart:
    with open("denoise_config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,

results = [
        f"denoise/{config['pool']}_error_R1.rds",
        f"denoise/{config['pool']}_asvs.fasta",
        f"denoise/{config['pool']}_asv_counts.tsv",
        f"denoise/multiqc_reports/{config['pool']}_denoise_metrics_mqc.out",
    ]
if is_paired():
   results.append(f"denoise/{config['pool']}_error_R2.rds")

rule all:
    input:
        results

rule dada2_learn_errors:
    """
    learn pool-level error profiles
    """
    input:
        manifest=config["manifest"],
    output:
        error=f"denoise/{config['pool']}_error_R{{dir}}.rds",
        error_fig=f"denoise/plots/{config['pool']}_dada2_error_estimation_R{{dir}}.png",
    log:
        e=f"{LOG_PREFIX}/dada2_learn_errors_R{{dir}}.e",
        o=f"{LOG_PREFIX}/dada2_learn_errors_R{{dir}}.o",
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 16
    resources:
        mem_mb=32 * 1024,
    script:
        "../scripts/denoise/dada2_learn_errors.R"

def get_infer_asv_fastq(wildcards):
    tmp = MANIFEST.loc[wildcards.sample, f"R{wildcards.dir}"]
    return tmp

rule dada2_infer_asvs:
    input:
        fastq=get_infer_asv_fastq,
        error=f"denoise/{config['pool']}_error_R{{dir}}.rds",
    output:
        derep=temp("denoise/dada2/{sample}_derep_R{dir}.rds"),
        dada=temp("denoise/dada2/{sample}_dada_R{dir}.rds"),
    params:
        pooling = "none"
    log:
        e=f"{LOG_PREFIX}/dada2_infer_asvs_{{sample}}_R{{dir}}.e",
        o=f"{LOG_PREFIX}/dada2_infer_asvs_{{sample}}_R{{dir}}.o",
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 16
    resources:
        mem_mb=32 * 1024,
    script:
        "../scripts/denoise/dada2_infer_asvs.R"


rule dada2_count_asvs:
    input:
        get_inputs_for_asv_counting
    output:
        merged=temp("denoise/dada2/{sample}_merged.rds"),
        seqtab=temp("denoise/dada2/{sample}_asv_seqtab.rds"),
    log:
        e=f"{LOG_PREFIX}/dada2_count_asvs_{{sample}}.e",
        o=f"{LOG_PREFIX}/dada2_count_asvs_{{sample}}.o",
    params:
        is_paired = is_paired()
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 1
    resources:
        mem_mb=8 * 1024,
    script:
        "../scripts/denoise/dada2_count_asvs.R"

use rule dada2_count_asvs as dada2_count_asvs_se with:
    output:
        seqtab=temp("denoise/dada2/{sample}_asv_seqtab-se.rds"),


rule dada2_remove_chimeras:
    input:
        seqtab=expand("denoise/dada2/{sample}_asv_seqtab{se}.rds", sample=SAMPLES, se="" if is_paired() else "-se"),
    output:
        seqtab=f"denoise/{config['pool']}_asv_seqtab.tsv",
        counts=f"denoise/{config['pool']}_asv_counts.tsv",
        fasta=f"denoise/{config['pool']}_asvs.fasta",
    log:
        e=f"{LOG_PREFIX}/dada2_remove_chimeras.e",
        o=f"{LOG_PREFIX}/dada2_remove_chimeras.o",
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 16
    resources:
        mem_mb=32 * 1024,
    script:
        "../scripts/denoise/dada2_remove_chimeras.R"


rule collect_dada2_sample_metrics:
    input:
        get_inputs_for_asv_counting,
        merged="denoise/dada2/{sample}_merged.rds",
    output:
        metrics=temp("denoise/dada2/{sample}_dada2_metrics.tsv"),
    params:
        is_paired = is_paired()
    log:
        e=f"{LOG_PREFIX}/collect_dada2_sample_metrics_{{sample}}.e",
        o=f"{LOG_PREFIX}/collect_dada2_sample_metrics_{{sample}}.o",
    resources:
        mem_mb=4 * 1024,
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    script:
        "../scripts/denoise/collect_dada2_sample_metrics.R"

use rule collect_dada2_sample_metrics as collect_dada2_sample_metrics_se with:
    input:
        get_inputs_for_asv_counting
    output:
        metrics=temp("denoise/dada2/{sample}_dada2_metrics-se.tsv"),


rule aggregate_metrics:
    input:
        sample_metrics=expand(
            "denoise/dada2/{sample}_dada2_metrics{se}.tsv",
            sample=SAMPLES,
            se="" if is_paired() else "-se"
        ),
        seq_counts=f"denoise/{config['pool']}_asv_seqtab.tsv",
    output:
        metrics=f"denoise/{config['pool']}_denoise_metrics.tsv",
    log:
        e=f"{LOG_PREFIX}/aggregate_metrics.e",
        o=f"{LOG_PREFIX}/aggregate_metrics.o",
    resources:
        mem_mb=4 * 1024,
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    script:
        "../scripts/denoise/aggregate_metrics.R"


rule make_denoise_metrics_multiqc_report:
    input:
        metrics=f"denoise/{config['pool']}_denoise_metrics.tsv",
    output:
        report=f"denoise/multiqc_reports/{config['pool']}_denoise_metrics_mqc.out",
    threads: 1
    log:
        e=f"{LOG_PREFIX}/make_denoise_metrics_multiqc_report.e",
        o=f"{LOG_PREFIX}/make_denoise_metrics_multiqc_report.o",
    resources:
        mem_mb=4 * 1024,
    shell:
        """
        echo -e "# id: 'dada2_metrics'" >  {output.report}
        echo -e "# plot_type: 'table'" >>  {output.report}
        echo -e "# section_name: 'DADA2 Denoising Metrics'" >>  {output.report}
        cat {input.metrics} >> {output.report}
        """
