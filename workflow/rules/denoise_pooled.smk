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

assert config["pooling"] != "none", "pooled workflow selected but config specifies pooling=none"

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
    threads: 8
    resources:
        mem_mb=8 * 1024,
    script:
        "../scripts/denoise/dada2_learn_errors.R"


rule dada2_infer_pooled_asvs:
    input:
        manifest=config["manifest"],
        error=f"denoise/{config['pool']}_error_R{{dir}}.rds",
    output:
        derep=f"denoise/dada2/{config['pool']}_derep_R{{dir}}.rds",
        dada=f"denoise/dada2/{config['pool']}_dada_R{{dir}}.rds",
    log:
        e=f"{LOG_PREFIX}/dada2_infer_asvs_{config['pool']}_R{{dir}}.e",
        o=f"{LOG_PREFIX}/dada2_infer_asvs_{config['pool']}_R{{dir}}.o",
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    params:
        pooling = config["pooling"]
    threads: 8
    resources:
        mem_mb=8 * 1024,
    script:
        "../scripts/denoise/dada2_infer_asvs.R"


rule dada2_postprocess:
    input:
        unpack(get_inputs_for_asv_counting)
    output:
        seqtab=f"denoise/{config['pool']}_asv_seqtab.tsv",
        counts=f"denoise/{config['pool']}_asv_counts.tsv",
        fasta=f"denoise/{config['pool']}_asvs.fasta",
        metrics=f"denoise/{config['pool']}_denoise_metrics.tsv",
        #merged=f"denoise/dada2/{config['pool']}_merged.rds",
    log:
        e=f"{LOG_PREFIX}/dada2_postprocess_{config['pool']}.e",
        o=f"{LOG_PREFIX}/dada2_postprocess_{config['pool']}.o",
    params:
        is_paired = is_paired()
    container:
        "docker://ghcr.io/vdblab/dada2:1.20.0"
    threads: 1
    resources:
        mem_mb=8 * 1024,
    script:
        "../scripts/denoise/dada2_pooled_postprocess.R"


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
