import os
import math
import traceback
import numpy as np
from contextlib import redirect_stderr


include: "common.smk"


envvars:
    "TMPDIR",


min_version("6.0")

validate(config, os.path.join(workflow.basedir, "../config/preprocess.schema.yaml"))

MANIFEST = load_manifest(config["manifest"], None)
SAMPLES = get_samples_from_manifest(MANIFEST)

LOG_PREFIX = "logs/preprocess"


onstart:
    with open("preprocess_config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,

def get_num_shards_for_sample(sample):
    """Calculate number of shards based on file size (1-20 shards)"""
    try:
        fq1_size = os.path.getsize(MANIFEST.loc[sample, "R1"])
        total_size = fq1_size
        if is_paired():
            fq2_size = os.path.getsize(MANIFEST.loc[sample, "R2"])
            total_size += fq2_size
        
        size_gb = total_size / (1024**3)
        return max(1, min(20, math.ceil(size_gb))) #20 seemed like a reasonable upper limit
    except:
        return 1

def get_fastq_list(wildcards):
    fqs = [MANIFEST.loc[wildcards.sample, "R1"]]
    if is_paired():
        fqs.append(MANIFEST.loc[wildcards.sample, "R2"])
    else:
        assert np.isnan(
            MANIFEST.loc[wildcards.sample, "R2"]
        ), f"Manifest contains R2 for sample {wildcards.sample}, but config lib_layout is not 'paired'"
    return fqs


def get_fastqs_to_trim(wildcards):
    """return a list of fastqs for trimming depending on whether paired library and
    whether primers have already been removed
    """
    if config["removeprimers"]:
        if is_paired():
            return [
                "preprocess/primers-removed/{sample}_noprimers_R1.fastq.gz",
                "preprocess/primers-removed/{sample}_noprimers_R2.fastq.gz",
            ]
        else:
            return ["preprocess/primers-removed-se/{sample}_noprimers_R1.fastq.gz"]
    else:
        # we keep the same basename so that fastqc works well for both
        if is_paired():
            return [
                "tmp/{sample}_noprimers_R1.fastq.gz",
                "tmp/{sample}_noprimers_R2.fastq.gz",
            ]
        else:
            return ["tmp/{sample}_noprimers_R1.fastq.gz"]


#        return get_fastq_list(wildcards)


def get_all_outputs(wildcards):
    if is_paired():
        results = [
            f"preprocess/{config['pool']}_manifest.tsv",
            f"preprocess/{config['pool']}_missing_or_incomplete.tsv",
            f"preprocess/{config['pool']}_failed_samples.tsv",
            f"preprocess/multiqc_reports/{config['pool']}_dada2_trimming_report_mqc.out",
            f"preprocess/multiqc_reports/{config['pool']}_adapter_contam_report_mqc.out",
        ]
        if config["removeprimers"]:
            results.extend(
                expand(
                    "preprocess/primers-removed/{sample}_noprimers_R{d}.fastq.gz",
                    sample=SAMPLES,
                    d=[1, 2],
                )
            )
    else:
        results = [
            f"preprocess/{config['pool']}_manifest.tsv",
            f"preprocess/{config['pool']}_missing_or_incomplete.tsv",
            f"preprocess/{config['pool']}_failed_samples.tsv",
            f"preprocess/multiqc_reports/{config['pool']}_dada2_trimming_report_mqc.out",
            f"preprocess/multiqc_reports/{config['pool']}_adapter_contam_report_mqc.out",
        ]
        if config["removeprimers"]:
            results.extend(
                expand(
                    "preprocess/primers-removed-se/{sample}_noprimers_R1.fastq.gz",
                    sample=SAMPLES,
                )
            )
    return results


rule all:
    input:
        get_all_outputs,


rule cp_if_not_removing_primers:
    """ make a local copy of input data that doesn't require primer removal.
    FastQC doesn't allow custom output file names, so in order to keep the names consistent we copy them here.
    symlinking isn't reliable on our cluster.
    I attempted to add a series of renaming commands after fastqc was run on the raw files, but
    because the outputs are compressed, the contents still contain the raw (unreliable) file names.
    This is a little bit of overhead but simplifies things greatly downstream.

    NOTE: if you (for some reason) rename the outputs, make sure to redefine get_fastqs_to_trim()!
    """
    input:
        get_fastq_list,
    output:
        R1=temp("tmp/{sample}_noprimers_R1.fastq.gz"),
        R2=temp("tmp/{sample}_noprimers_R2.fastq.gz"),
    shell:
        """
    # we know the inputs are ordered
    counter=1
    touch {output.R2}
    for f in {input}
    do
    cp $f tmp/{wildcards.sample}_noprimers_R${{counter}}.fastq.gz
    counter=$((counter + 1))
    done
   """

checkpoint split_fastq:
    input:
        reads=get_fastq_list,
    output:
        directory("preprocess/chunks/{sample}")
    run:
        import subprocess
        num_shards = get_num_shards_for_sample(wildcards.sample)
        
        # If only 1 shard, just symlink the original files
        if num_shards == 1:
            os.makedirs(output[0], exist_ok=True)
            shell(f"ln -sf $(realpath {input.reads[0]}) {output[0]}/chunk_00_R1.fastq.gz")
            if len(input.reads) > 1:
                shell(f"ln -sf $(realpath {input.reads[1]}) {output[0]}/chunk_00_R2.fastq.gz")
        else:
            # Split the files
            os.makedirs(output[0], exist_ok=True)
            
            # Count reads for splitting
            if input.reads[0].endswith('.gz'):
                read_count_cmd = f"zcat {input.reads[0]} | wc -l"
            else:
                read_count_cmd = f"wc -l < {input.reads[0]}"
            
            total_lines = int(subprocess.check_output(read_count_cmd, shell=True).decode().strip())
            reads_per_shard = max(1, math.ceil((total_lines // 4) / num_shards))
            lines_per_shard = reads_per_shard * 4
            
            # Split R1
            if input.reads[0].endswith('.gz'):
                shell(f"zcat {input.reads[0]} | split -l {lines_per_shard} -d --additional-suffix='_R1.fastq' - {output[0]}/chunk_")
            else:
                shell(f"split -l {lines_per_shard} -d --additional-suffix='_R1.fastq' {input.reads[0]} {output[0]}/chunk_")
            
            # Split R2 if paired-end
            if len(input.reads) > 1:
                if input.reads[1].endswith('.gz'):
                    shell(f"zcat {input.reads[1]} | split -l {lines_per_shard} -d --additional-suffix='_R2.fastq' - {output[0]}/chunk_")
                else:
                    shell(f"split -l {lines_per_shard} -d --additional-suffix='_R2.fastq' {input.reads[1]} {output[0]}/chunk_")
            
            # Gzip the chunks
            shell(f"gzip {output[0]}/chunk_*.fastq")

rule remove_primers_chunk:
    input:
        R1="preprocess/chunks/{sample}/chunk_{chunk}_R1.fastq.gz",
        R2="preprocess/chunks/{sample}/chunk_{chunk}_R2.fastq.gz"
    output:
        R1=temp("preprocess/chunks-processed/{sample}/chunk_{chunk}_R1.fastq.gz"),
        R2=temp("preprocess/chunks-processed/{sample}/chunk_{chunk}_R2.fastq.gz")
    params:
        outstring=lambda wc, input, output: (
            " -o {} -p {} ".format(output.R1, output.R2)
            if os.path.exists(str(input.R2))
            else " -o {}".format(output.R1)
        ),
        primerstring=lambda wc, input, output: (
            " -g correct={F} -g flipped={R} -G flipped={F} -G correct={R}".format(
                F=config["primer_F"], R=config["primer_R"]
            )
            if os.path.exists(str(input.R2))
            else " -g correct={F} -g flipped={R} ".format(
                F=config["primer_F"], R=config["primer_R"]
            )
        ),
        input_files=lambda wc, input: f"{input.R1} {input.R2}" if os.path.exists(str(input.R2)) else str(input.R1),
        preprocess_min_read_length=config["preprocess_min_read_length"],
    container:
        "docker://ghcr.io/vdblab/cutadapt:3.7"
    threads: 4
    resources:
        mem_mb=4 * 1024,
    log:
        o=f"{LOG_PREFIX}/cutadapt_remove_primers_{{sample}}_chunk_{{chunk}}.o",
        e=f"{LOG_PREFIX}/cutadapt_remove_primers_{{sample}}_chunk_{{chunk}}.e",
    shell:
        """
        cutadapt \
            -Z \
            --cores {threads} \
            {params.primerstring} \
            -e 1.5 \
            --minimum-length {params.preprocess_min_read_length} \
            {params.outstring} \
            --rename='{{header}}  {{adapter_name}}' \
            {params.input_files} \
            > {log.o} 2>> {log.e}
        """

def aggregate_primer_chunks(wildcards):
    checkpoint_output = checkpoints.split_fastq.get(sample=wildcards.sample).output[0]
    chunks = glob_wildcards(os.path.join(checkpoint_output, "chunk_{chunk}_R1.fastq.gz")).chunk
    
    result = {
        "R1": expand("preprocess/chunks-processed/{sample}/chunk_{chunk}_R1.fastq.gz", 
                    sample=wildcards.sample, chunk=chunks)
    }
    
    # Check if R2 exists for first chunk to determine if paired-end
    if os.path.exists(f"{checkpoint_output}/chunk_{chunks[0]}_R2.fastq.gz"):
        result["R2"] = expand("preprocess/chunks-processed/{sample}/chunk_{chunk}_R2.fastq.gz", 
                            sample=wildcards.sample, chunk=chunks)
    
    return result

rule remove_primers:
    """
    Concatenate processed chunks back together
    """
    input:
        unpack(aggregate_primer_chunks)
    output:
        R1=temp("preprocess/primers-removed/{sample}_noprimers_R1.fastq.gz"),
        R2=temp("preprocess/primers-removed/{sample}_noprimers_R2.fastq.gz"),
    shell:
        """
        # Concatenate R1 files
        cat {input.R1} > {output.R1}
        
        # Concatenate R2 files if they exist
        if [ ! -z "{input.R2}" ]; then
            cat {input.R2} > {output.R2}
        else
            touch {output.R2}
        fi
        """

use rule remove_primers as remove_primers_se with:
    output:
        R1=temp("preprocess/primers-removed-se/{sample}_noprimers_R1.fastq.gz"),


rule quality_trim:
    input:
        reads=get_fastqs_to_trim,
    output:
        R1="preprocess/trimmed/{sample}_trimmed_R1.fastq.gz",
        R2="preprocess/trimmed/{sample}_trimmed_R2.fastq.gz",
        figpath="preprocess/plots/{sample}_quality_profile_mqc.png",
        stats=temp("preprocess/dada2/{sample}_quality_trimming_stats.txt"),
    params:
        trunclen_R1=config["trunclen_R1"],
        trunclen_R2=config["trunclen_R2"],
        trimLeft=config["dada2_trimLeft"],
        trimRight=config["dada2_trimRight"],
        maxLen=config["dada2_maxLen"],
        maxN=config["dada2_maxN"],
        maxEE=config["dada2_maxEE"],
        truncQ=config["dada2_truncQ"],
        rmphix=config["dada2_rm.phix"],
        rmlowcomplex=config["dada2_rm.lowcomplex"],
        figdir_name="dada2/figures/",
        min_asv_len=config["min_asv_len"],
        is_paired=is_paired(),
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


use rule quality_trim as quality_trim_se with:
    output:
        R1="preprocess/trimmed-se/{sample}_trimmed_R1.fastq.gz",
        figpath="preprocess/plots/{sample}_quality_profile-se_mqc.png",
        stats=temp("preprocess/dada2/{sample}_quality_trimming_stats-se.txt"),


rule output_manifest:
    """
    Output a file with sample ids extracted from the oligo files,
    and trimmed R1 and R1 fastq files
    """
    input:
        expand(
            "preprocess/trimmed{se}/{sample}_trimmed_R{dir}.fastq.gz",
            se="" if is_paired() else "-se",
            sample=SAMPLES,
            dir=[1, 2] if is_paired() else [1],
        ),
    output:
        manifest=f"preprocess/{config['pool']}_manifest.tsv",
        missing=f"preprocess/{config['pool']}_missing_or_incomplete.tsv",
    log:
        e=f"{LOG_PREFIX}/output_manifest.e",
    run:
        with open(log.e, "w") as ef, redirect_stderr(ef):
            try:
                fastq_template = os.path.join(
                    os.getcwd(),
                    os.path.dirname(input[0]),
                    "{sample}_trimmed_R{dir}.fastq.gz",
                )
                write_manifest_and_missing(
                    SAMPLES,
                    fastq_template,
                    output["manifest"],
                    output["missing"],
                    paired=is_paired(),
                )
            except Exception as e:
                traceback.print_exc(file=ef)


rule sample_fastqc_report:
    input:
        get_fastqs_to_trim,
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
            {input} \
            2> {log.e} > {log.o}
        """


use rule sample_fastqc_report as sample_fastqc_report_se with:
    output:
        report_R1="preprocess/fastqc_reports-se/{sample}_noprimers_R1_fastqc.html",
        zip_R1="preprocess/fastqc_reports-se/{sample}_noprimers_R1_fastqc.zip",


def get_fastqc_outputs(wildcards):
    if is_paired():
        return [
            "preprocess/fastqc_reports/{sample}_noprimers_R1_fastqc.zip",
            "preprocess/fastqc_reports/{sample}_noprimers_R2_fastqc.zip",
        ]
    else:
        return ["preprocess/fastqc_reports-se/{sample}_noprimers_R1_fastqc.zip"]


rule get_contam_from_fastqc:
    """
    this prints the overrepresented sequences section with awk and sed,
    then removes all the "No Hit" entries that are likely biologically overrepresented,
    such as the amplicon's primer binding site, non-varying regions between V4-V5, etc
    """
    input:
        get_fastqc_outputs,
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


def get_trim_reports(wildcards):
    if is_paired():
        return expand(
            "preprocess/dada2/{sample}_quality_trimming_stats.txt", sample=SAMPLES
        )
    else:
        return expand(
            "preprocess/dada2/{sample}_quality_trimming_stats-se.txt", sample=SAMPLES
        )


rule make_dada2_trimming_multiqc_report:
    """
    This is its own rule purely so we can mark the input files as temporary,
    which will make generating unit tests possible
    """
    input:
        get_trim_reports,
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
        cat {input} | grep -v "reads.in" >> {output.report} 2>> {log.e}
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
